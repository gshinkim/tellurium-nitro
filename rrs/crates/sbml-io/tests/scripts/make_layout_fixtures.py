"""Writes the Layout test files in tests/data with python-libsbml:
a Level 3 Version 1 model with a layout that uses every kind of glyph, and
the same model in Level 2 Version 4, where the layout lives in the model's
annotation.

    python make_layout_fixtures.py <tests/data folder>
"""
import sys, os, libsbml

out = sys.argv[1]
def build(level, version):
    ns = libsbml.LayoutPkgNamespaces(level, version, 1)
    doc = libsbml.SBMLDocument(ns)
    if level == 3:
        doc.setPackageRequired("layout", False)
    m = doc.createModel()
    m.setId("glycolysis_fragment")
    c = m.createCompartment(); c.setId("cell"); c.setSize(1); c.setConstant(True); c.setSpatialDimensions(3)
    for sid, amt in [("Glc", 10), ("G6P", 0), ("ATP", 5), ("ADP", 0), ("HK", 1)]:
        s = m.createSpecies(); s.setId(sid); s.setCompartment("cell"); s.setInitialConcentration(amt)
        s.setHasOnlySubstanceUnits(False); s.setBoundaryCondition(False); s.setConstant(False)
    p = m.createParameter(); p.setId("k1"); p.setValue(0.1); p.setConstant(True)
    r = m.createReaction(); r.setId("J0"); r.setReversible(False); r.setFast(False)
    for sid, ref in [("Glc", "sr_glc"), ("ATP", "sr_atp")]:
        sr = r.createReactant(); sr.setSpecies(sid); sr.setId(ref); sr.setStoichiometry(1); sr.setConstant(True)
    for sid, ref in [("G6P", "sr_g6p"), ("ADP", "sr_adp")]:
        sr = r.createProduct(); sr.setSpecies(sid); sr.setId(ref); sr.setStoichiometry(1); sr.setConstant(True)
    mr = r.createModifier(); mr.setSpecies("HK"); mr.setId("sr_hk")
    kl = r.createKineticLaw(); kl.setMath(libsbml.parseL3Formula("k1*HK*Glc*ATP"))

    lp = m.getPlugin("layout")
    lay = lp.createLayout(); lay.setId("layout1"); lay.setDimensions(libsbml.Dimensions(ns, 400, 300))
    cg = lay.createCompartmentGlyph(); cg.setId("cg_cell"); cg.setCompartmentId("cell")
    cg.setBoundingBox(libsbml.BoundingBox(ns, "bb_cell", 5, 5, 390, 290))
    pos = {"Glc": (20, 20), "G6P": (300, 20), "ATP": (20, 200), "ADP": (300, 200), "HK": (160, 120)}
    for sid, (x, y) in pos.items():
        g = lay.createSpeciesGlyph(); g.setId("sg_" + sid); g.setSpeciesId(sid)
        g.setBoundingBox(libsbml.BoundingBox(ns, "bb_" + sid, x, y, 60, 30))
        t = lay.createTextGlyph(); t.setId("tg_" + sid); t.setOriginOfTextId(sid); t.setGraphicalObjectId("sg_" + sid)
        t.setBoundingBox(libsbml.BoundingBox(ns, "bbt_" + sid, x, y, 60, 30))
    rg = lay.createReactionGlyph(); rg.setId("rg_J0"); rg.setReactionId("J0")
    seg = rg.getCurve().createLineSegment()
    seg.setStart(libsbml.Point(ns, 180, 100)); seg.setEnd(libsbml.Point(ns, 200, 100))
    roles = {"Glc": ("sr_glc", libsbml.SPECIES_ROLE_SUBSTRATE), "ATP": ("sr_atp", libsbml.SPECIES_ROLE_SIDESUBSTRATE),
             "G6P": ("sr_g6p", libsbml.SPECIES_ROLE_PRODUCT), "ADP": ("sr_adp", libsbml.SPECIES_ROLE_SIDEPRODUCT),
             "HK": ("sr_hk", libsbml.SPECIES_ROLE_ACTIVATOR)}
    for sid, (ref, role) in roles.items():
        srg = rg.createSpeciesReferenceGlyph(); srg.setId("srg_" + sid); srg.setSpeciesGlyphId("sg_" + sid)
        srg.setSpeciesReferenceId(ref); srg.setRole(role)
        x, y = pos[sid]
        if sid in ("Glc", "G6P"):
            b = srg.getCurve().createCubicBezier()
            b.setStart(libsbml.Point(ns, 190, 100)); b.setEnd(libsbml.Point(ns, x + 30, y + 15))
            b.setBasePoint1(libsbml.Point(ns, 190, 60.5)); b.setBasePoint2(libsbml.Point(ns, x + 30, 60.25))
        else:
            l = srg.getCurve().createLineSegment()
            l.setStart(libsbml.Point(ns, 190, 100)); l.setEnd(libsbml.Point(ns, x + 30, y + 15))
    go = lay.createAdditionalGraphicalObject(); go.setId("legend")
    go.setBoundingBox(libsbml.BoundingBox(ns, "bb_legend", 10, 260, 80, 20))
    gg = lay.createGeneralGlyph(); gg.setId("gg_k1"); gg.setReferenceId("k1")
    gg.setBoundingBox(libsbml.BoundingBox(ns, "bb_k1", 160, 60, 40, 20))
    ref = gg.createReferenceGlyph(); ref.setId("ref_k1"); ref.setGlyphId("rg_J0"); ref.setRole("parameter")
    l = ref.getCurve().createLineSegment(); l.setStart(libsbml.Point(ns, 180, 70)); l.setEnd(libsbml.Point(ns, 190, 100))
    sub = libsbml.TextGlyph(ns); sub.setId("tg_k1"); sub.setText("k1")
    sub.setBoundingBox(libsbml.BoundingBox(ns, "bb_tg_k1", 160, 60, 40, 20))
    gg.addSubGlyph(sub)
    return doc

def check(d, name):
    # libsbml's consistency check objects to a curve and a bounding box on
    # the same glyph, which libsbml itself writes; check the file reads back.
    d = libsbml.readSBMLFromString(libsbml.writeSBMLToString(d))
    errs = [d.getError(i).getMessage() for i in range(d.getNumErrors())
            if d.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR]
    assert not errs, (name, errs[:3])

for level, version in [(3, 1), (2, 4)]:
    doc = build(level, version)
    name = f"layout-l{level}v{version}"
    check(doc, name)
    libsbml.writeSBMLToFile(doc, os.path.join(out, name + ".xml"))
print("wrote", out)
