from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions

def _config_draw_options(draw_options: dict, size: tuple[int]) -> tuple:
    '''
    Configure drawer and drawing options given provided options.
    '''
    drawer = Draw.MolDraw2DSVG(*size)
    _draw_options = drawer.drawOptions()
    for k, v in draw_options.items():
        if not hasattr(_draw_options, k):
            raise ValueError(f"Select from {dir(_draw_options)}")
        elif callable(getattr(_draw_options, k)) and v is not None:
            getattr(_draw_options, k)(v)
        elif callable(getattr(_draw_options, k)) and v is None:
            getattr(_draw_options, k)()
        else:
            setattr(_draw_options, k, v)
    return drawer, _draw_options

def draw_reaction(rxn: str, sub_img_size: tuple = (300, 200), use_smiles: bool = True, draw_options: dict = {}) -> str:
    '''
    Draw reaction to svg string

    Args
    -----
    rxn:str
        Reaction SMARTS
    sub_img_size:tuple
        width by height
    use_smiles:bool
        If True, is more explicit about double
        bond location in drawing
    draw_options:dict
        Key-value pairs to set fields or call functions of rdkit.Chem.Draw.rdMolDraw2D.drawOptions.
        For functions that take no arguments, set the value to None. For functions that take arguments,
        set the value to the argument. For fields, set the value to the desired value.

        Examples:
            draw_options = {
                'addAtomIndices': True,
                'useBWAtomPalette': None,
                'setBackgroundColour': (1.0, 1.0, 1.0),
                'setHighlightColour': (1.0, 0.0, 0.0),
                'highlightBondWidthMultiplier': 2.0,
            }

        Below docstring of options cf:
        https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
            - property addAtomIndices
                adds atom indices to drawings. Default False.

            - property addBondIndices
                adds bond indices to drawings. Default False.

            - property addStereoAnnotation
                adds R/S and E/Z to drawings. Default False.

            - property additionalAtomLabelPadding
                additional padding to leave around atom labels. Expressed as a fraction of the font size.

            - property annotationFontScale
                Scale of font for atom and bond annotation relative to atomlabel font. Default=0.75.

            - property atomHighlightsAreCircles
            forces atom highlights always to be circles.Default (false) is to put ellipses roundlonger labels.

            - property atomLabelDeuteriumTritium
                labels deuterium as D and tritium as T

            - property atomLabels
                maps indices to atom labels

            - property atomRegions
                regions to outline

            - property baseFontSize
                relative size of font. Defaults to 0.6. -1 means use default.

            - property bondLineWidth
                if positive, this overrides the default line width for bonds

            - property bracketsAroundAtomLists
                Whether to put brackets round atom lists in query atoms. Default is true.

            - property centreMoleculesBeforeDrawing
                Moves the centre of the drawn molecule to (0,0).Default False.

            - property circleAtoms
                
            - property clearBackground
                clear the background before drawing a molecule

            - property comicMode
                simulate hand-drawn lines for bonds. When combined with a font like Comic-Sans or Comic-Neue, this gives xkcd-like drawings. Default is false.

            - property continuousHighlight
                
            - property drawMolsSameScale
                when drawing multiple molecules with DrawMolecules, forces them to use the same scale. Default is true.

            - property drawingExtentsInclude
                Drawing extents are computed taking into account only selected DrawElement items. Default=DrawElement.ALL

            - property dummiesAreAttachments

            - property dummyIsotopeLabels
                adds isotope labels on dummy atoms. Default True.

            - property explicitMethyl
                Draw terminal methyls explictly. Default is false.

            - property fillHighlights

            - property fixedBondLength
                If > 0.0, fixes bond length to this number of pixelsunless that would make it too big. Default -1.0 meansno fix. If both set, fixedScale takes precedence.

            - property fixedFontSize
                font size in pixels. default=-1 means not fixed. If set, always used irrespective of scale, minFontSize and maxFontSize.

            - property fixedScale
                If > 0.0, fixes scale to that fraction of width ofdraw window. Default -1.0 means adjust scale to fit.

            - property flagCloseContactsDist

            - property fontFile
                Font file for use with FreeType text drawer. Can also be BuiltinTelexRegular (the default) or BuiltinRobotoRegular.

            - getAnnotationColour((MolDrawOptions)self) → object :
                method returning the annotation colour

                C++ signature :
                boost::python::api::object getAnnotationColour(RDKit::MolDrawOptions)

            - getAtomNoteColour((MolDrawOptions)self) → object :
                method returning the atom note colour

                C++ signature :
                boost::python::api::object getAtomNoteColour(RDKit::MolDrawOptions)

            - getBackgroundColour((MolDrawOptions)self) → object :
                method returning the background colour

                C++ signature :
                boost::python::api::object getBackgroundColour(RDKit::MolDrawOptions)

            - getBondNoteColour((MolDrawOptions)self) → object :
                method returning the bond note colour

                C++ signature :
                boost::python::api::object getBondNoteColour(RDKit::MolDrawOptions)

            - getHighlightColour((MolDrawOptions)self) → object :
                method returning the highlight colour

                C++ signature :
                boost::python::api::object getHighlightColour(RDKit::MolDrawOptions)

            - getLegendColour((MolDrawOptions)self) → object :
                method returning the legend colour

                C++ signature :
                boost::python::api::object getLegendColour(RDKit::MolDrawOptions)

            - getQueryColour((MolDrawOptions)self) → object :
                method returning the query colour

                C++ signature :
                boost::python::api::object getQueryColour(RDKit::MolDrawOptions)

            - getSymbolColour((MolDrawOptions)self) → object :
                method returning the symbol colour

                C++ signature :
                boost::python::api::object getSymbolColour(RDKit::MolDrawOptions)

            - getVariableAttachmentColour((MolDrawOptions)self) → object :
                method for getting the colour of variable attachment points

                C++ signature :
                boost::python::api::object getVariableAttachmentColour(RDKit::MolDrawOptions)

            - property highlightBondWidthMultiplier
                What to multiply default bond width by for highlighting bonds. Default-8.

            - property highlightRadius
                Default radius for highlight circles.

            - property includeAtomTags
                include atom tags in output

            - property includeChiralFlagLabel
                add a molecule annotation with “ABS” if the chiral flag is set. Default is false.

            - property includeMetadata
                When possible, include metadata about molecules and reactions to allow them to be reconstructed. Default is true.

            - property includeRadicals
                include radicals in the drawing (it can be useful to turn this off for reactions and queries). Default is true.

            - property isotopeLabels
                adds isotope labels on non-dummy atoms. Default True.

            - property legendFontSize
                font size in pixels of the legend (if drawn)

            - property legendFraction
                fraction of the draw panel to be used for the legend if present

            - property maxFontSize
                maximum font size in pixels. default=40, -1 means no maximum.

            - property minFontSize
                minimum font size in pixels. default=6, -1 means no minimum.

            - property multiColourHighlightStyle
                Either 'CircleAndLine' or 'Lasso', to control style ofmulti-coloured highlighting in DrawMoleculeWithHighlights.Default is CircleAndLine.

            - property multipleBondOffset
                offset for the extra lines in a multiple bond as a fraction of mean bond length

            - property noAtomLabels
                disables inclusion of atom labels in the rendering

            - property padding
                Fraction of empty space to leave around molecule. Default=0.05.

            - property prepareMolsBeforeDrawing
                call prepareMolForDrawing() on each molecule passed to DrawMolecules()

            - property reagentPadding
                Fraction of empty space to leave around each component of a reaction drawing. Default=0.0.

            - property rotate
                Rotates molecule about centre by this number of degrees,

            - property scaleBondWidth
                Scales the width of drawn bonds using image scaling.

            - property scaleHighlightBondWidth
                Scales the width of drawn highlighted bonds using image scaling.

            - property scalingFactor
                scaling factor for pixels->angstrom when auto scalingbeing used. Default is 20.

            - setAnnotationColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the annotation colour

                C++ signature :
                void setAnnotationColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setAtomNoteColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the atom note colour

                C++ signature :
                void setAtomNoteColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setAtomPalette((MolDrawOptions)self, (AtomPairsParameters)cmap) → None :
                sets the palette for atoms and bonds from a dictionary mapping ints to 3-tuples

                C++ signature :
                void setAtomPalette(RDKit::MolDrawOptions {lvalue},boost::python::api::object)

            - setBackgroundColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the background colour

                C++ signature :
                void setBackgroundColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setBondNoteColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the bond note colour

                C++ signature :
                void setBondNoteColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setHighlightColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the highlight colour

                C++ signature :
                void setHighlightColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setLegendColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the legend colour

                C++ signature :
                void setLegendColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setQueryColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the query colour

                C++ signature :
                void setQueryColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setSymbolColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the symbol colour

                C++ signature :
                void setSymbolColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setVariableAttachmentColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the colour of variable attachment points

                C++ signature :
                void setVariableAttachmentColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - property showAllCIPCodes
                show all defined CIP codes (no hiding!). Default False.

            - property simplifiedStereoGroupLabel
                if all specified stereocenters are in a single StereoGroup, show a molecule-level annotation instead of the individual labels. Default is false.

            - property singleColourWedgeBonds
                if true wedged and dashed bonds are drawn using symbolColour rather than inheriting their colour from the atoms. Default is false.

            - property splitBonds
            - property standardColoursForHighlightedAtoms
                If true, highlighted hetero atoms are drawn in standard colours rather than black. Default=False

            - property unspecifiedStereoIsUnknown
                if true, double bonds with unspecified stereo are drawn crossed, potential stereocenters with unspecified stereo are drawn with a wavy bond. Default is false.

            - updateAtomPalette((MolDrawOptions)self, (AtomPairsParameters)cmap) → None :
                updates the palette for atoms and bonds from a dictionary mapping ints to 3-tuples

                C++ signature :
                void updateAtomPalette(RDKit::MolDrawOptions {lvalue},boost::python::api::object)

            - useAvalonAtomPalette((MolDrawOptions)self) → None :
                use the Avalon renderer palette for atoms and bonds

                C++ signature :
                void useAvalonAtomPalette(RDKit::MolDrawOptions {lvalue})

            - useBWAtomPalette((MolDrawOptions)self) → None :
                use a black and white palette for atoms and bonds

                C++ signature :
                void useBWAtomPalette(RDKit::MolDrawOptions {lvalue})

            - useCDKAtomPalette((MolDrawOptions)self) → None :
                use the CDK palette for atoms and bonds

                C++ signature :
                void useCDKAtomPalette(RDKit::MolDrawOptions {lvalue})

            - property useComplexQueryAtomSymbols
                replace any atom, any hetero, any halo queries with complex query symbols A, Q, X, M, optionally followed by H if hydrogen is included (except for AH, which stays *). Default is true

            - useDefaultAtomPalette((MolDrawOptions)self) → None :
                use the default colour palette for atoms and bonds

                C++ signature :
                void useDefaultAtomPalette(RDKit::MolDrawOptions {lvalue})

            - property useMolBlockWedging
                If the molecule came from a MolBlock, prefer the wedging information that provides. If false, use RDKit rules. Default false

            - property variableAtomRadius
                radius value to use for atoms involved in variable attachment points.

            - property variableBondWidthMultiplier
                what to multiply standard bond width by for variable attachment points.
    
    Returns
    -------
    :str
        SVG string of drawn reaction
    '''
    _, _draw_options = _config_draw_options(draw_options, sub_img_size)
    rxn = rdChemReactions.ReactionFromSmarts(rxn, useSmiles=use_smiles)
    return Draw.ReactionToImage(rxn, useSVG=True, subImgSize=sub_img_size, drawOptions=_draw_options)

def draw_molecule(molecule: str | Chem.Mol, size: tuple = (200, 200), highlight_atoms: tuple = tuple(), draw_options: dict = {}, legend: str = '') -> str:
    '''
    Draw molecule to svg string

    Args
    ----
    mol:str | Chem.Mol
        Molecule
    size:tuple
        (width, height)
    highlight_atoms:tuple
        Atom indices to highlight
    draw_options:dict
        Key-value pairs to set fields or call functions of rdkit.Chem.Draw.rdMolDraw2D.drawOptions. For
        functions that take no arguments, set the value to None. For functions that take arguments, set the value to the argument.
        For fields, set the value to the desired value.

        Examples:
            draw_options = {
                'addAtomIndices': True,
                'useBWAtomPalette': None,
                'setBackgroundColour': (1.0, 1.0, 1.0),
                'setHighlightColour': (1.0, 0.0, 0.0),
                'highlightBondWidthMultiplier': 2.0,
            }

        Below docstring of options cf:
        https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
            - property addAtomIndices
                adds atom indices to drawings. Default False.

            - property addBondIndices
                adds bond indices to drawings. Default False.

            - property addStereoAnnotation
                adds R/S and E/Z to drawings. Default False.

            - property additionalAtomLabelPadding
                additional padding to leave around atom labels. Expressed as a fraction of the font size.

            - property annotationFontScale
                Scale of font for atom and bond annotation relative to atomlabel font. Default=0.75.

            - property atomHighlightsAreCircles
            forces atom highlights always to be circles.Default (false) is to put ellipses roundlonger labels.

            - property atomLabelDeuteriumTritium
                labels deuterium as D and tritium as T

            - property atomLabels
                maps indices to atom labels

            - property atomRegions
                regions to outline

            - property baseFontSize
                relative size of font. Defaults to 0.6. -1 means use default.

            - property bondLineWidth
                if positive, this overrides the default line width for bonds

            - property bracketsAroundAtomLists
                Whether to put brackets round atom lists in query atoms. Default is true.

            - property centreMoleculesBeforeDrawing
                Moves the centre of the drawn molecule to (0,0).Default False.

            - property circleAtoms

            - property clearBackground
                clear the background before drawing a molecule

            - property comicMode
                simulate hand-drawn lines for bonds. When combined with a font like Comic-Sans or Comic-Neue, this gives xkcd-like drawings. Default is false.

            - property continuousHighlight

            - property drawMolsSameScale
                when drawing multiple molecules with DrawMolecules, forces them to use the same scale. Default is true.

            - property drawingExtentsInclude
                Drawing extents are computed taking into account only selected DrawElement items. Default=DrawElement.ALL

            - property dummiesAreAttachments

            - property dummyIsotopeLabels
                adds isotope labels on dummy atoms. Default True.

            - property explicitMethyl
                Draw terminal methyls explictly. Default is false.

            - property fillHighlights

            - property fixedBondLength
                If > 0.0, fixes bond length to this number of pixelsunless that would make it too big. Default -1.0 meansno fix. If both set, fixedScale takes precedence.

            - property fixedFontSize
                font size in pixels. default=-1 means not fixed. If set, always used irrespective of scale, minFontSize and maxFontSize.

            - property fixedScale
                If > 0.0, fixes scale to that fraction of width ofdraw window. Default -1.0 means adjust scale to fit.

            - property flagCloseContactsDist

            - property fontFile
                Font file for use with FreeType text drawer. Can also be BuiltinTelexRegular (the default) or BuiltinRobotoRegular.

            - getAnnotationColour((MolDrawOptions)self) → object :
                method returning the annotation colour

                C++ signature :
                boost::python::api::object getAnnotationColour(RDKit::MolDrawOptions)

            - getAtomNoteColour((MolDrawOptions)self) → object :
                method returning the atom note colour

                C++ signature :
                boost::python::api::object getAtomNoteColour(RDKit::MolDrawOptions)

            - getBackgroundColour((MolDrawOptions)self) → object :
                method returning the background colour

                C++ signature :
                boost::python::api::object getBackgroundColour(RDKit::MolDrawOptions)

            - getBondNoteColour((MolDrawOptions)self) → object :
                method returning the bond note colour

                C++ signature :
                boost::python::api::object getBondNoteColour(RDKit::MolDrawOptions)

            - getHighlightColour((MolDrawOptions)self) → object :
                method returning the highlight colour

                C++ signature :
                boost::python::api::object getHighlightColour(RDKit::MolDrawOptions)

            - getLegendColour((MolDrawOptions)self) → object :
                method returning the legend colour

                C++ signature :
                boost::python::api::object getLegendColour(RDKit::MolDrawOptions)

            - getQueryColour((MolDrawOptions)self) → object :
                method returning the query colour

                C++ signature :
                boost::python::api::object getQueryColour(RDKit::MolDrawOptions)

            - getSymbolColour((MolDrawOptions)self) → object :
                method returning the symbol colour

                C++ signature :
                boost::python::api::object getSymbolColour(RDKit::MolDrawOptions)

            - getVariableAttachmentColour((MolDrawOptions)self) → object :
                method for getting the colour of variable attachment points

                C++ signature :
                boost::python::api::object getVariableAttachmentColour(RDKit::MolDrawOptions)

            - property highlightBondWidthMultiplier
                What to multiply default bond width by for highlighting bonds. Default-8.

            - property highlightRadius
                Default radius for highlight circles.

            - property includeAtomTags
                include atom tags in output

            - property includeChiralFlagLabel
                add a molecule annotation with “ABS” if the chiral flag is set. Default is false.

            - property includeMetadata
                When possible, include metadata about molecules and reactions to allow them to be reconstructed. Default is true.

            - property includeRadicals
                include radicals in the drawing (it can be useful to turn this off for reactions and queries). Default is true.

            - property isotopeLabels
                adds isotope labels on non-dummy atoms. Default True.

            - property legendFontSize
                font size in pixels of the legend (if drawn)

            - property legendFraction
                fraction of the draw panel to be used for the legend if present

            - property maxFontSize
                maximum font size in pixels. default=40, -1 means no maximum.

            - property minFontSize
                minimum font size in pixels. default=6, -1 means no minimum.

            - property multiColourHighlightStyle
                Either 'CircleAndLine' or 'Lasso', to control style ofmulti-coloured highlighting in DrawMoleculeWithHighlights.Default is CircleAndLine.

            - property multipleBondOffset
                offset for the extra lines in a multiple bond as a fraction of mean bond length

            - property noAtomLabels
                disables inclusion of atom labels in the rendering

            - property padding
                Fraction of empty space to leave around molecule. Default=0.05.

            - property prepareMolsBeforeDrawing
                call prepareMolForDrawing() on each molecule passed to DrawMolecules()

            - property reagentPadding
                Fraction of empty space to leave around each component of a reaction drawing. Default=0.0.

            - property rotate
                Rotates molecule about centre by this number of degrees,

            - property scaleBondWidth
                Scales the width of drawn bonds using image scaling.

            - property scaleHighlightBondWidth
                Scales the width of drawn highlighted bonds using image scaling.

            - property scalingFactor
                scaling factor for pixels->angstrom when auto scalingbeing used. Default is 20.

            - setAnnotationColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the annotation colour

                C++ signature :
                void setAnnotationColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setAtomNoteColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the atom note colour

                C++ signature :
                void setAtomNoteColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setAtomPalette((MolDrawOptions)self, (AtomPairsParameters)cmap) → None :
                sets the palette for atoms and bonds from a dictionary mapping ints to 3-tuples

                C++ signature :
                void setAtomPalette(RDKit::MolDrawOptions {lvalue},boost::python::api::object)

            - setBackgroundColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the background colour

                C++ signature :
                void setBackgroundColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setBondNoteColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the bond note colour

                C++ signature :
                void setBondNoteColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setHighlightColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the highlight colour

                C++ signature :
                void setHighlightColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setLegendColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the legend colour

                C++ signature :
                void setLegendColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setQueryColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the query colour

                C++ signature :
                void setQueryColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setSymbolColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the symbol colour

                C++ signature :
                void setSymbolColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - setVariableAttachmentColour((MolDrawOptions)self, (tuple)tpl) → None :
                method for setting the colour of variable attachment points

                C++ signature :
                void setVariableAttachmentColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)

            - property showAllCIPCodes
                show all defined CIP codes (no hiding!). Default False.

            - property simplifiedStereoGroupLabel
                if all specified stereocenters are in a single StereoGroup, show a molecule-level annotation instead of the individual labels. Default is false.

            - property singleColourWedgeBonds
                if true wedged and dashed bonds are drawn using symbolColour rather than inheriting their colour from the atoms. Default is false.

            - property splitBonds
            - property standardColoursForHighlightedAtoms
                If true, highlighted hetero atoms are drawn in standard colours rather than black. Default=False

            - property unspecifiedStereoIsUnknown
                if true, double bonds with unspecified stereo are drawn crossed, potential stereocenters with unspecified stereo are drawn with a wavy bond. Default is false.

            - updateAtomPalette((MolDrawOptions)self, (AtomPairsParameters)cmap) → None :
                updates the palette for atoms and bonds from a dictionary mapping ints to 3-tuples

                C++ signature :
                void updateAtomPalette(RDKit::MolDrawOptions {lvalue},boost::python::api::object)

            - useAvalonAtomPalette((MolDrawOptions)self) → None :
                use the Avalon renderer palette for atoms and bonds

                C++ signature :
                void useAvalonAtomPalette(RDKit::MolDrawOptions {lvalue})

            - useBWAtomPalette((MolDrawOptions)self) → None :
                use a black and white palette for atoms and bonds

                C++ signature :
                void useBWAtomPalette(RDKit::MolDrawOptions {lvalue})

            - useCDKAtomPalette((MolDrawOptions)self) → None :
                use the CDK palette for atoms and bonds

                C++ signature :
                void useCDKAtomPalette(RDKit::MolDrawOptions {lvalue})

            - property useComplexQueryAtomSymbols
                replace any atom, any hetero, any halo queries with complex query symbols A, Q, X, M, optionally followed by H if hydrogen is included (except for AH, which stays *). Default is true

            - useDefaultAtomPalette((MolDrawOptions)self) → None :
                use the default colour palette for atoms and bonds

                C++ signature :
                void useDefaultAtomPalette(RDKit::MolDrawOptions {lvalue})

            - property useMolBlockWedging
                If the molecule came from a MolBlock, prefer the wedging information that provides. If false, use RDKit rules. Default false

            - property variableAtomRadius
                radius value to use for atoms involved in variable attachment points.

            - property variableBondWidthMultiplier
                what to multiply standard bond width by for variable attachment points.
    
    Returns
    -------
    :str
        SVG string of drawn molecule
    '''
    if type(molecule) is str:
        mol = Chem.MolFromSmiles(molecule)

        # Catch failed MolFromSmiles
        if mol is None: 
            mol = Chem.MolFromSmiles(molecule, sanitize=False)
    elif type(molecule) is Chem.Mol:
        mol = molecule

    drawer, _ = _config_draw_options(draw_options, size)
    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, legend=legend)
    drawer.FinishDrawing()
    img = drawer.GetDrawingText()

    return img

if __name__ == '__main__':
    draw_molecule('CCO', draw_options={'addAtomIndices': True})
    draw_reaction('[C:1]=[O:2].[N:3]>>[C:1][N:3].[O:2]', draw_options={'addAtomIndices': True})