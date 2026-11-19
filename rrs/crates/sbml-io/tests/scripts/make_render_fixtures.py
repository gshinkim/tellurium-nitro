"""Writes the Render test files in tests/data with python-libsbml: a Level 3
Version 1 model with a layout, a local render information on the layout and a
global one, between them using every part of Render (colours, both kinds of
gradient, line endings, styles by role, type and id, and every primitive);
and the same document converted by libsbml to Level 2 Version 4, where layout
and Render live in annotations.

    python make_render_fixtures.py <tests/data folder>
"""
import sys, os, libsbml

out = sys.argv[1]
RA = libsbml.RelAbsVector


def model_and_layout(level, version):
    ns = libsbml.LayoutPkgNamespaces(level, version, 1)
    doc = libsbml.SBMLDocument(ns)
    if level == 3:
        doc.setPackageRequired("layout", False)
        doc.enablePackage(libsbml.RenderExtension.getXmlnsL3V1V1(), "render", True)
        doc.setPackageRequired("render", False)
    else:
        doc.enablePackage(libsbml.RenderExtension.getXmlnsL2(), "render", True)
    m = doc.createModel()
    m.setId("render_fixture")
    c = m.createCompartment(); c.setId("cell"); c.setSize(1); c.setConstant(True); c.setSpatialDimensions(3)
    for sid in ["S", "P", "E"]:
        s = m.createSpecies(); s.setId(sid); s.setCompartment("cell"); s.setInitialConcentration(1)
        s.setHasOnlySubstanceUnits(False); s.setBoundaryCondition(False); s.setConstant(False)
    r = m.createReaction(); r.setId("J0"); r.setReversible(False); r.setFast(False)
    sr = r.createReactant(); sr.setSpecies("S"); sr.setId("sr_s"); sr.setStoichiometry(1); sr.setConstant(True)
    sr = r.createProduct(); sr.setSpecies("P"); sr.setId("sr_p"); sr.setStoichiometry(1); sr.setConstant(True)
    mr = r.createModifier(); mr.setSpecies("E"); mr.setId("sr_e")

    lay = m.getPlugin("layout").createLayout()
    lay.setId("layout1"); lay.setDimensions(libsbml.Dimensions(ns, 400, 300))
    cg = lay.createCompartmentGlyph(); cg.setId("cg_cell"); cg.setCompartmentId("cell")
    cg.setBoundingBox(libsbml.BoundingBox(ns, "bb_cell", 5, 5, 390, 290))
    pos = {"S": (20, 120), "P": (300, 120), "E": (160, 20)}
    for sid, (x, y) in pos.items():
        g = lay.createSpeciesGlyph(); g.setId("sg_" + sid); g.setSpeciesId(sid)
        g.setBoundingBox(libsbml.BoundingBox(ns, "bb_" + sid, x, y, 60, 30))
        t = lay.createTextGlyph(); t.setId("tg_" + sid); t.setOriginOfTextId(sid); t.setGraphicalObjectId("sg_" + sid)
        t.setBoundingBox(libsbml.BoundingBox(ns, "bbt_" + sid, x, y, 60, 30))
    rg = lay.createReactionGlyph(); rg.setId("rg_J0"); rg.setReactionId("J0")
    rg.setBoundingBox(libsbml.BoundingBox(ns, "bb_J0", 186, 131, 8, 8))
    for sid, ref, role in [("S", "sr_s", libsbml.SPECIES_ROLE_SUBSTRATE), ("P", "sr_p", libsbml.SPECIES_ROLE_PRODUCT),
                           ("E", "sr_e", libsbml.SPECIES_ROLE_MODIFIER)]:
        srg = rg.createSpeciesReferenceGlyph(); srg.setId("srg_" + sid); srg.setSpeciesGlyphId("sg_" + sid)
        srg.setSpeciesReferenceId(ref); srg.setRole(role)
        x, y = pos[sid]
        l = srg.getCurve().createLineSegment()
        l.setStart(libsbml.Point(ns, 190, 135)); l.setEnd(libsbml.Point(ns, x + 30, y + 15))
    t = lay.createTextGlyph(); t.setId("tg_note"); t.setText("A note")
    t.setBoundingBox(libsbml.BoundingBox(ns, "bb_note", 10, 260, 120, 20))
    return doc, lay


def render_info(ri, ns_render):
    """Fills a render information with one of everything."""
    ri.setProgramName("make_render_fixtures.py")
    ri.setProgramVersion("1")
    ri.setBackgroundColor("#F8F9FAFF")
    for cid, value in [("black", "#000000"), ("fill", "#EEF6FF"), ("teal", "#007777CC"), ("none_like", "#FFFFFF00")]:
        c = ri.createColorDefinition(); c.setId(cid); c.setColorValue(value)
    lg = ri.createLinearGradientDefinition(); lg.setId("lin"); lg.setSpreadMethod("reflect")
    lg.setPoint1(RA(0, 0), RA(0, 0)); lg.setPoint2(RA(0, 100), RA(0, 100))
    for off, col in [(0, "#FFFFFF"), (100, "fill")]:
        s = lg.createGradientStop(); s.setOffset(RA(0, off)); s.setStopColor(col)
    rgd = ri.createRadialGradientDefinition(); rgd.setId("rad")
    rgd.setCenter(RA(0, 50), RA(0, 50)); rgd.setRadius(RA(0, 50)); rgd.setFocalPoint(RA(0, 40), RA(0, 40))
    s = rgd.createGradientStop(); s.setOffset(RA(0, 0)); s.setStopColor("#FF0000")
    s = rgd.createGradientStop(); s.setOffset(RA(0, 100)); s.setStopColor("black")

    le = ri.createLineEnding(); le.setId("arrow"); le.setEnableRotationalMapping(True)
    le.setBoundingBox(libsbml.BoundingBox(ns_render, "bb_arrow", -12, -5, 12, 10))
    p = le.getGroup().createPolygon(); p.setFill("black")
    for x, y in [(0, 0), (100, 50), (0, 100)]:
        p.createPoint().setCoordinates(RA(0, x), RA(0, y))
    le = ri.createLineEnding(); le.setId("dot"); le.setEnableRotationalMapping(False)
    le.setBoundingBox(libsbml.BoundingBox(ns_render, "bb_dot", -4, -4, 8, 8))
    e = le.getGroup().createEllipse(); e.setCenter2D(RA(0, 50), RA(0, 50)); e.setRadii(RA(0, 50), RA(0, 50))
    e.setFill("teal")
    return ri


def styles(ri, local):
    st = ri.createStyle("species")
    st.addType("SPECIESGLYPH")
    g = st.getGroup()
    g.setStroke("black"); g.setStrokeWidth(1.5); g.setFillColor("fill")
    g.setFontFamily("sans-serif"); g.setFontSize(RA(12, 0)); g.setFontWeight(libsbml.FONT_WEIGHT_BOLD)
    g.setFontStyle(libsbml.FONT_STYLE_ITALIC); g.setTextAnchor(libsbml.H_TEXTANCHOR_MIDDLE)
    g.setVTextAnchor(libsbml.V_TEXTANCHOR_MIDDLE)
    rect = g.createRectangle()
    rect.setCoordinatesAndSize(RA(0, 0), RA(0, 0), RA(0, 0), RA(0, 100), RA(0, 100))
    rect.setRadiusX(RA(8, 0)); rect.setRadiusY(RA(8, 0))
    rect.addDash(5); rect.addDash(4)

    st = ri.createStyle("reaction_lines")
    st.addRole("product")
    st.addRole("substrate")
    g = st.getGroup(); g.setStroke("black"); g.setStrokeWidth(2); g.setEndHead("arrow")

    st = ri.createStyle("modifier_lines")
    st.addRole("modifier")
    g = st.getGroup(); g.setStroke("teal"); g.setStrokeWidth(1.5); g.setEndHead("dot")
    g.addDash(6); g.addDash(4)

    st = ri.createStyle("compartments")
    st.addType("COMPARTMENTGLYPH")
    g = st.getGroup(); g.setStroke("black"); g.setFillColor("rad")
    e = g.createEllipse(); e.setCenter2D(RA(0, 50), RA(0, 50)); e.setRadii(RA(0, 50), RA(-2, 50))

    if local:
        st = ri.createStyle("just_P")
        st.addId("sg_P")
        g = st.getGroup(); g.setStroke("#C80000"); g.setFillColor("lin")
        p = g.createPolygon()
        for x, y in [(0, 50), (25, 0), (75, 0), (100, 50), (75, 100), (25, 100)]:
            p.createPoint().setCoordinates(RA(0, x), RA(0, y))
        cb = p.createCubicBezier()
        cb.setCoordinates(RA(0, 0), RA(0, 50))
        cb.setBasePoint1(RA(0, 10), RA(0, 90))
        cb.setBasePoint2(RA(0, 5), RA(0, 70))

        st = ri.createStyle("just_E")
        st.addId("sg_E")
        st.addId("tg_E")
        g = st.getGroup(); g.setStroke("black")
        inner = g.createGroup(); inner.setStroke("teal"); inner.setStrokeWidth(3)
        c = inner.createCurve()
        c.createPoint().setCoordinates(RA(0, 0), RA(0, 0))
        c.createPoint().setCoordinates(RA(0, 100), RA(0, 100))
        c.setEndHead("arrow")
        t = g.createText(); t.setText("E"); t.setCoordinates(RA(0, 50), RA(0, 50))
        t.setFontSize(RA(10, 0)); t.setTextAnchor(libsbml.H_TEXTANCHOR_MIDDLE)
        im = g.createImage(); im.setImageReference("enzyme.png")
        im.setCoordinates(RA(0, 0), RA(0, 0)); im.setWidth(RA(0, 100)); im.setHeight(RA(0, 100))

        st = ri.createStyle("note")
        st.addType("TEXTGLYPH")
        g = st.getGroup(); g.setStroke("#1A1A1A"); g.setFontFamily("serif"); g.setFontSize(RA(14, 0))
        g.setTextAnchor(libsbml.H_TEXTANCHOR_START)


def build(level, version, name):
    doc, lay = model_and_layout(level, version)
    rns = libsbml.LayoutPkgNamespaces(level, version, 1)
    lp = doc.getModel().getPlugin("layout")
    gri = lp.getListOfLayouts().getPlugin("render").createGlobalRenderInformation()
    gri.setId("global1"); gri.setName("Global style")
    render_info(gri, rns)
    styles(gri, local=False)
    lri = lay.getPlugin("render").createLocalRenderInformation()
    lri.setId("local1"); lri.setReferenceRenderInformationId("global1")
    render_info(lri, rns)
    styles(lri, local=True)
    path = os.path.join(out, name)
    libsbml.writeSBMLToFile(doc, path)
    # Checked as read back from the file.
    back = libsbml.readSBMLFromFile(path)
    back.checkConsistency()
    log = back.getErrorLog()
    for i in range(log.getNumErrors()):
        e = log.getError(i)
        if e.getSeverity() >= libsbml.LIBSBML_SEV_ERROR:
            print(name, "ERROR", e.getErrorId(), e.getMessage()[:200])


build(3, 1, "render-l3v1.xml")
build(2, 4, "render-l2v4.xml")
print("written")
