"""Accept all tracked changes in a .docx by rewriting the XML directly.

Written because the skill's accept_changes.py reports success on a LibreOffice
timeout while copying the input through unchanged, and headless soffice macro
dispatch silently no-ops in this container.

Comments (comments.xml and their anchors) are deliberately PRESERVED -- only
revision markup is resolved.
"""
import re, shutil, sys, zipfile
from pathlib import Path
from lxml import etree

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
def w(t): return f"{{{W}}}{t}"

# Revision containers whose CONTENT survives (the change is kept).
UNWRAP = {w("ins"), w("moveTo")}
# Revision containers whose content goes away entirely.
DROP = {w("del"), w("moveFrom")}
# Formatting-revision records: drop the record, keep current formatting.
DROP_MARKER = {
    w("pPrChange"), w("rPrChange"), w("tblPrChange"), w("tcPrChange"),
    w("trPrChange"), w("sectPrChange"), w("tblGridChange"),
    w("numberingChange"), w("cellIns"), w("cellDel"), w("cellMerge"),
    w("moveFromRangeStart"), w("moveFromRangeEnd"),
    w("moveToRangeStart"), w("moveToRangeEnd"),
}

# Anchors that must SURVIVE a dropped deletion: removing them orphans a
# reviewer comment or a cross-reference target.
RESCUE_BARE = {w("commentRangeStart"), w("commentRangeEnd"),
               w("bookmarkStart"), w("bookmarkEnd")}
RESCUE_WRAP = {w("commentReference")}


def rescue_from(el):
    """Pull comment/bookmark anchors out of a subtree about to be deleted.

    Returns them in document order, ready to be spliced in where the
    subtree stood. A commentReference must sit inside a run, so it gets a
    fresh empty one.
    """
    saved = []
    for n in el.iter():
        if n.tag in RESCUE_BARE:
            saved.append(n)
        elif n.tag in RESCUE_WRAP:
            r = etree.Element(w("r"))
            r.append(n)
            saved.append(r)
    for n in saved:
        # Detach from the doomed subtree without disturbing sibling order.
        if n.getparent() is not None and n.getparent() is not el:
            pass
    return saved


def para_mark_deleted(p):
    """True if this paragraph's MARK was deleted (=> merge into the next one)."""
    pPr = p.find(w("pPr"))
    if pPr is None:
        return False
    rPr = pPr.find(w("rPr"))
    return rPr is not None and rPr.find(w("del")) is not None

def accept_tree(root):
    stats = {"ins": 0, "del": 0, "moveTo": 0, "moveFrom": 0, "fmt": 0, "joins": 0}

    # 1. Note which paragraphs lose their mark BEFORE we strip the w:del records.
    marked = [p for p in root.iter(w("p")) if para_mark_deleted(p)]

    # 2. Drop deleted content and formatting-revision records.
    for el in list(root.iter()):
        if el.tag in DROP:
            # A w:del inside pPr/rPr is the paragraph-mark record, handled above.
            parent = el.getparent()
            if parent is not None and parent.tag == w("rPr") \
               and parent.getparent() is not None and parent.getparent().tag == w("pPr"):
                parent.remove(el)
                continue
            stats["del" if el.tag == w("del") else "moveFrom"] += 1
            saved = rescue_from(el)
            gp = el.getparent()
            idx = list(gp).index(el)
            gp.remove(el)
            for k, node in enumerate(saved):
                gp.insert(idx + k, node)
        elif el.tag in DROP_MARKER:
            stats["fmt"] += 1
            el.getparent().remove(el)

    # 3. Unwrap kept insertions, preserving document order.
    for el in list(root.iter()):
        if el.tag in UNWRAP:
            stats["ins" if el.tag == w("ins") else "moveTo"] += 1
            parent = el.getparent()
            idx = list(parent).index(el)
            for child in reversed(list(el)):
                parent.insert(idx, child)
            parent.remove(el)

    # 4. Merge each mark-deleted paragraph into the following one.
    #    Reverse order so chains (p1->p2->p3) collapse correctly.
    for p in reversed(marked):
        parent = p.getparent()
        if parent is None:
            continue
        sibs = list(parent)
        i = sibs.index(p)
        nxt = next((s for s in sibs[i + 1:] if s.tag == w("p")), None)
        if nxt is None or sibs[i + 1] is not nxt:
            continue  # nothing to merge into, or a table/etc. intervenes
        content = [c for c in p if c.tag != w("pPr")]
        anchor = nxt.find(w("pPr"))
        pos = (list(nxt).index(anchor) + 1) if anchor is not None else 0
        for c in content:
            nxt.insert(pos, c)
            pos += 1
        parent.remove(p)
        stats["joins"] += 1

    return stats

def main(src, dst):
    src, dst = Path(src), Path(dst)
    shutil.copy2(src, dst)
    zin = zipfile.ZipFile(src)
    targets = [n for n in zin.namelist()
               if re.match(r"word/(document|header\d*|footer\d*|footnotes|endnotes)\.xml$", n)]
    total = {}
    parts = {}
    for name in targets:
        root = etree.fromstring(zin.read(name))
        s = accept_tree(root)
        for k, v in s.items():
            total[k] = total.get(k, 0) + v
        parts[name] = etree.tostring(root, xml_declaration=True,
                                     encoding="UTF-8", standalone=True)
    # Rewrite the archive, replacing only the parts we touched.
    with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            if item.filename.endswith("/"):
                continue
            data = parts.get(item.filename, zin.read(item.filename))
            zout.writestr(item, data)
    zin.close()
    print(f"{src.name} -> {dst.name}")
    for k, v in total.items():
        print(f"  {k:9s} {v}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
