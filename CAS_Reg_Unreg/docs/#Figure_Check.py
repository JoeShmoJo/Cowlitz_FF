#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure provenance: keep the memo, the staging folder, and the script outputs
from drifting apart.

THE PROBLEM THIS EXISTS TO SOLVE
--------------------------------
Every memo figure exists in three places. The script writes it to
<project>/output/, a copy is staged in docs/figures_combined/, and a third copy
is embedded in the .docx. Nothing links the three. A script can be re-run and
the document keeps the old picture, and there is no way to see that by looking
at the document. It has already happened once: the Coweeman figure was
regenerated and the stale staging copy went into the memo anyway.

Three copies with no link is the hazard. This script removes the drift by
making the staging folder a BUILD PRODUCT rather than a hand-kept pile, and by
checking the document against the live output.

    MODE = "check"   report every figure in the memo against its live source
    MODE = "stage"   rebuild figures_combined/ from the live sources
    MODE = "both"    stage, then check

STAGED FILES CARRY THE FIGURE NUMBER
    Fig_5-8__unreg_reg_final_uncertainty.png

That is the real protection. You cannot paste the wrong figure when the file
name says which figure it is, and you cannot paste a stale one because staging
copies from the live output every time it runs.

WHAT "STALE" MEANS HERE
    The embedded image is compared to the live source by bytes first. Bytes
    differing is not proof of staleness, because matplotlib does not write
    byte-identical files across versions or dpi settings. So a byte mismatch
    falls through to a pixel comparison when Pillow is available. Same pixels
    means the figure is current and was merely re-rendered. Different pixels
    means the document is behind the script. Without Pillow the check reports
    BYTES-DIFFER and says it could not confirm, rather than crying wolf.

RUN IT
    python3 "#Figure_Check.py"

WHEN A FIGURE DOES NOT MATCH
    The check does not guess which side is newer. File timestamps look like the
    obvious signal and they are not, because a fresh git clone stamps every
    file with the checkout time. Instead the embedded image is hunted for
    across the staging folder and both output trees, and the report names the
    file it actually came from. DIFFERS with a named file usually means a stale
    staging copy was pasted in over the live one, which is the failure this
    whole script exists to catch. NO ORIGIN means the picture in the document
    is reproducible from nothing at all, which is worse.

Nothing here writes to the .docx. Checking is read only, staging only touches
docs/figures_combined/.
"""

import os
import re
import csv
import sys
import shutil
import hashlib
import zipfile
import datetime

os.chdir(os.path.dirname(os.path.abspath(__file__)))

MODE = "both"

DOCX = "MEMO_CAS_Combined_FlowFrequency_2026_03_09.docx"
MANIFEST = "figure_manifest.csv"
STAGE_DIR = "figures_combined"

# Remove staged files that the manifest does not name. Off by default, since
# the folder still holds a few hand made images.
STAGE_PRUNE = False

# Mean absolute pixel difference, 0 to 255, below which two renderings of the
# same figure are treated as the same picture. Word downscales images on
# insert, so a same-content pair lands around 1 to 3 while a genuinely
# different plot lands above 10. Measured on this memo: Figure 7-1, which Word
# shrank from 2025x870 to 1428x613, and the Appendix C figures, which really
# are out of date.
PIXEL_TOLERANCE = 6.0

# Where to look for the true origin of a figure that does not match its
# manifest source. Relative to this folder.
SEARCH_ROOTS = ["figures_combined",
                "../output",
                "../../CAS_Unreg_FF/output"]

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def md5_bytes(b):
    return hashlib.md5(b).hexdigest()


def md5_file(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def mtime(path):
    t = datetime.datetime.fromtimestamp(os.path.getmtime(path))
    return t.strftime("%Y-%m-%d %H:%M")


def read_manifest(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if not r.get("figure"):
                continue
            rows.append({"figure": r["figure"].strip(),
                         "source": (r.get("source") or "").strip(),
                         "script": (r.get("script") or "").strip(),
                         "note": (r.get("note") or "").strip()})
    return rows


def docx_figures(docx_path):
    """Pair each embedded image with the Figure caption that follows it.

    Word stores the picture and its caption as separate paragraphs, so the
    caption is found by walking paragraphs in document order and attaching the
    pending image to the next paragraph whose text starts with "Figure".
    """
    z = zipfile.ZipFile(docx_path)
    rels = dict(re.findall(r'Id="([^"]+)"[^>]*Target="media/([^"]+)"',
                           z.read("word/_rels/document.xml.rels").decode("utf-8")))
    xml = z.read("word/document.xml").decode("utf-8")
    body = xml[xml.index("<w:body"):]
    out = []
    pending = []
    for para in re.findall(r"<w:p[ >].*?</w:p>", body, re.S):
        for rid in re.findall(r'r:embed="([^"]+)"', para):
            if rid in rels:
                pending.append(rels[rid])
        text = "".join(re.findall(r"<w:t(?:\s[^>]*)?>(.*?)</w:t>", para, re.S))
        text = re.sub(r"\s+", " ", text).strip()
        m = re.match(r"^(Figure [A-Za-z0-9-]+)", text)
        if m and pending:
            media = pending[0]
            out.append({"figure": m.group(1),
                        "media": media,
                        "bytes": z.read("word/media/" + media)})
            pending = []
    z.close()
    return out


def same_pixels(blob, path):
    """Compare an embedded image to a file on disk.

    Returns "same", "rescaled", "different", or None when it cannot be decided.

    Word re-compresses and downscales pictures on insert, so an embedded copy is
    routinely a smaller rendering of the same figure. Comparing raw bytes calls
    that stale, and a checker that cries wolf gets ignored, which is how the
    stale figure got into the memo in the first place. So a size mismatch at the
    same aspect ratio is resolved by scaling the larger image down to the
    smaller and comparing the mean absolute pixel difference. A re-render of the
    same plot lands near zero. A genuinely different plot lands far above the
    threshold, since axis limits, lines and text all move.
    """
    try:
        from PIL import Image
        import io as _io
    except Exception:
        return None
    try:
        a = Image.open(_io.BytesIO(blob)).convert("RGB")
        b = Image.open(path).convert("RGB")
    except Exception:
        return None
    if a.size == b.size:
        return "same" if a.tobytes() == b.tobytes() else "different"
    ar = (a.size[0] / float(a.size[1])) / (b.size[0] / float(b.size[1]))
    if not (0.98 < ar < 1.02):
        return "different"
    target = a.size if a.size[0] * a.size[1] < b.size[0] * b.size[1] else b.size
    a = a.resize(target, Image.LANCZOS)
    b = b.resize(target, Image.LANCZOS)
    n = target[0] * target[1] * 3
    diff = sum(abs(x - y) for x, y in zip(a.tobytes(), b.tobytes())) / float(n)
    return "rescaled" if diff <= PIXEL_TOLERANCE else "different"


def locate(blob):
    """Find where the memo's copy of an image actually came from.

    Better than guessing direction from file timestamps, which a fresh git
    clone destroys. If the embedded image matches a file in the staging folder
    but not the live script output, the staging copy is stale and was pasted
    in over the real one. If it matches nothing at all, the figure has no
    reproducible origin. Either way the answer is concrete and checkable.
    """
    exact, close = [], []
    h = md5_bytes(blob)
    for root in SEARCH_ROOTS:
        for dirpath, _dirs, files in os.walk(root):
            for fn in files:
                if not fn.lower().endswith((".png", ".jpg", ".jpeg")):
                    continue
                full = os.path.normpath(os.path.join(dirpath, fn))
                try:
                    if md5_file(full) == h:
                        exact.append(full)
                        continue
                except Exception:
                    continue
                # Only fall back to the pixel test when nothing matched byte
                # for byte. The frequency curves by duration look nearly alike,
                # so a loose pixel match alone will name the wrong one.
                try:
                    if same_pixels(blob, full) in ("same", "rescaled"):
                        close.append(full)
                except Exception:
                    pass
    return exact if exact else close


def stage(rows):
    print("STAGING  ->  %s/" % STAGE_DIR)
    if not os.path.isdir(STAGE_DIR):
        os.makedirs(STAGE_DIR)
    written = []
    for r in rows:
        if r["source"] in ("", "-"):
            continue
        if not os.path.exists(r["source"]):
            print("   MISSING SOURCE  %-11s %s" % (r["figure"], r["source"]))
            continue
        num = r["figure"].replace("Figure ", "").replace(" ", "")
        name = "Fig_%s__%s" % (num, os.path.basename(r["source"]))
        dest = os.path.join(STAGE_DIR, name)
        old = md5_file(dest) if os.path.exists(dest) else None
        shutil.copyfile(r["source"], dest)
        written.append(name)
        flag = "     " if old == md5_file(dest) else "  NEW"
        print("  %s %-11s %-52s  (source written %s)"
              % (flag, r["figure"], name, mtime(r["source"])))
    with open(os.path.join(STAGE_DIR, "MANIFEST.txt"), "w") as f:
        f.write("Staged by #Figure_Check.py on %s\n"
                "Do not hand edit anything named Fig_*. Re-run this script "
                "instead.\n\n" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        for r in rows:
            if r["source"] in ("", "-"):
                f.write("%-11s  (supplied by hand)  %s\n" % (r["figure"], r["note"]))
                continue
            if not os.path.exists(r["source"]):
                continue
            num = r["figure"].replace("Figure ", "").replace(" ", "")
            name = "Fig_%s__%s" % (num, os.path.basename(r["source"]))
            f.write("%-11s  %s\n              from %s\n              md5 %s  written %s\n"
                    % (r["figure"], name, r["source"],
                       md5_file(r["source"]), mtime(r["source"])))
    if STAGE_PRUNE:
        for f in sorted(os.listdir(STAGE_DIR)):
            if f.startswith("Fig_") and f not in written:
                os.remove(os.path.join(STAGE_DIR, f))
                print("   PRUNED  %s" % f)
    print("   %d figures staged, MANIFEST.txt written\n" % len(written))


def check(rows):
    if not os.path.exists(DOCX):
        sys.exit("Cannot find %s. Set DOCX at the top of this script." % DOCX)
    print("CHECKING %s  against the live script outputs" % DOCX)
    print("   document written %s\n" % mtime(DOCX))
    by_fig = dict((r["figure"], r) for r in rows)
    found = docx_figures(DOCX)
    counts = {}
    print("  %-11s %-9s %s" % ("figure", "verdict", "detail"))
    print("  " + "-" * 96)
    for f in found:
        r = by_fig.get(f["figure"])
        if r is None:
            v, d = "UNLISTED", "not in %s, add it or the figure goes unchecked" % MANIFEST
        elif r["source"] in ("", "-"):
            v, d = "BY HAND", r["note"] or "no script source, nothing to check"
        elif not os.path.exists(r["source"]):
            v, d = "NO SOURCE", r["source"]
        elif md5_bytes(f["bytes"]) == md5_file(r["source"]):
            v, d = "CURRENT", "%s  (%s)" % (os.path.basename(r["source"]),
                                            mtime(r["source"]))
        else:
            px = same_pixels(f["bytes"], r["source"])
            if px == "same":
                v, d = "CURRENT", ("%s  re-rendered, pixels identical"
                                   % os.path.basename(r["source"]))
            elif px == "rescaled":
                v, d = "CURRENT", ("%s  same picture, Word rescaled it on insert"
                                   % os.path.basename(r["source"]))
            elif px == "different":
                where = locate(f["bytes"])
                live = os.path.normpath(r["source"])
                where = [w for w in where if w != live]
                if not where:
                    v = "NO ORIGIN"
                    d = ("the copy in the document matches no file on disk, so "
                         "nothing regenerates it. Re-run %s and compare."
                         % (r["script"] or "its script"))
                else:
                    v = "DIFFERS"
                    d = ("document holds %s, the live output is %s -- decide "
                         "which is right, then re-run %s"
                         % (where[0], os.path.basename(r["source"]),
                            r["script"] or "its script"))
            else:
                v, d = "BYTES DIFFER", ("could not decode to compare pixels, "
                                        "install Pillow to resolve  (%s)"
                                        % os.path.basename(r["source"]))
        counts[v] = counts.get(v, 0) + 1
        print("  %-11s %-9s %s" % (f["figure"], v, d))

    missing = [r["figure"] for r in rows
               if r["figure"] not in set(x["figure"] for x in found)]
    if missing:
        print("\n  in the manifest but not found in the document: %s"
              % ", ".join(missing))
    print("\n  " + "   ".join("%s %d" % (k, v) for k, v in sorted(counts.items())))
    bad = counts.get("DIFFERS", 0) + counts.get("NO ORIGIN", 0)
    if bad:
        print("\n  *** %d figure(s) in the document do not match the live "
              "script output ***" % bad)
    print()


def main():
    rows = read_manifest(MANIFEST)
    print("=" * 100)
    print("FIGURE PROVENANCE CHECK   %d figures in %s" % (len(rows), MANIFEST))
    print("=" * 100 + "\n")
    if MODE in ("stage", "both"):
        stage(rows)
    if MODE in ("check", "both"):
        check(rows)


main()
