def __find_param(var, input_cats, split_pars_dict, split_cat_pars):
    from P2VV.RooFitWrappers import RooObject
    ws = RooObject.ws()
    split_cats = split_pars_dict.get( ws[ var.GetName() ], set() )
    if not var.isFundamental() :
        split = False
        for splVar, splCats in split_pars_dict.iteritems() :
            if var.dependsOn(splVar) : split = True
            split_cats |= splCats
        if not split : split_cats = set()
    if not split_cats :
        return var
    else :
        catLabels = [(cat.GetName(), cat.getLabel()) for cat in input_cats if cat in split_cats]
        catsStr = ';'.join(lab[1] for lab in catLabels)
        if len(catLabels) > 1 : catsStr = '{' + catsStr + '}'
        from P2VV.Utilities.General import getSplitPar
        splitVar = getSplitPar(var.GetName(), catsStr, split_cat_pars)
        assert splitVar, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: parameter "%s" is set to be constrained for category "%s", but it is not found in PDF'\
               % (var.GetName(), catsStr)
        from ROOT import RooRealVar, RooCategory, RooFormulaVar
        from P2VV.RooFitWrappers import RealVar, Category, FormulaVar
        wrappers = { RooRealVar : RealVar, RooCategory : Category, RooFormulaVar : FormulaVar }
        assert type(splitVar) in wrappers\
               , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: wrapping of RooFit object "%s" of type "%s" not implemented (yet)'\
                 % ( splitVar.GetName(), type(splitVar) )
        return wrappers[ type(splitVar) ]( Name = splitVar.GetName() )

def __make_wrapper(orig_params, split_model, input_cats, split_params_dict, split_cat_pars):
    params = [ ]
    for var in orig_params:
        params.append(__find_param(var, input_cats, split_params_dict, split_cat_pars))
    from P2VV.RooFitWrappers import ResolutionModel
    return ResolutionModel(Name = split_model.GetName()), params

def __orig_params(orig_params, name):
    params = filter(lambda e: e[0] in name, orig_params.iteritems())
    assert(len(params) == 1)
    return params[0][1]

def replacement_model(prototype, split_model, input_cats, split_pars_dict,
                      split_cat_pars, orig_res_params):
    from P2VV.Parameterizations.TimeResolution import TimeResolution
    from P2VV.RooFitWrappers import AddModel
    from ROOT import RooAddModel
    if isinstance(split_model, RooAddModel):
        models = []
        for model in split_model.pdfList():
            model, params = __make_wrapper(__orig_params(orig_res_params, model.GetName()),
                                           model, input_cats, split_pars_dict, split_cat_pars)
            model._target_()._parameters = params
            model['ConditionalObservables'] = prototype['ConditionalObservables']
            models.append(model)
        fractions = [__find_param(f, input_cats, split_pars_dict, split_cat_pars) for f in __orig_params(orig_res_params, split_model.GetName())]
        return TimeResolution(Model = AddModel(split_model.GetName(), Models = models, Fractions = fractions))
    else:
        model, params = __make_wrapper(orig_res_params.values()[0], split_model, input_cats,
                                       split_pars_dict, split_cat_pars)
        return TimeResolution(Model = model, Parameters = params)
