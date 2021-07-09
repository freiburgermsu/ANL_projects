# import the COBRA model
import cobra
import optlang
bigg_model_path = '..\COBRA function scripts\e_coli_core metabolism from BiGG.json'
#bigg_model_path = 'e_coli_core metabolism from BiGG.json'
model = cobra.io.load_json_model(bigg_model_path)

# import the modelseedpy packages
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.fbapkg.simplethermopkg import SimpleThermoPkg
#from modelseedpy.fbapkg.fullthermopkg import FullThermoPkg
import sys
package_path = '..\ModelSEEDpy\modelseedpy\\fbapkg'
sys.path.insert(0, package_path)
from fullthermopkg import FullThermoPkg

# ------------------------ test the BaseFBA Package ---------------------------------------

full_thermo = FullThermoPkg(model = model)    

# import the modelseed database content
from modelseedpy.core.fbahelper import FBAHelper
modelseed_db_path = '..\..\..\Biofilm growth code\GSWL code\ModelSEEDDatabase'
ms_api = FBAHelper.get_modelseed_db_api(modelseed_db_path)
        
def test_init():
    
    # assert results of the model 
    assert type(full_thermo.model) is cobra.core.model.Model
    assert type(full_thermo.name) is str
    assert type(full_thermo.variable_types) is dict
    assert type(full_thermo.constraint_types) is dict
    assert 'logconc' in list(full_thermo.variables.keys())
    assert 'dgerr' in list(full_thermo.variables.keys())
    assert 'potc' in list(full_thermo.constraints.keys())

    # reinstantiate the model and clear the variables
    bigg_model_path = 'e_coli_core metabolism from BiGG.json'
    model = cobra.io.load_json_model(bigg_model_path)
    
    
def test_build_package():
    
    # execute the function
    full_thermo.build_package(full_thermo.parameters)

    # assert results from creating the variables and constraints
    added_parameters = ["default_max_conc", "default_min_conc", "default_max_error", "custom_concentrations",
                        "custom_deltaG_error", "compartment_potential", "temperature", "filter",'modelseed_path',
                        'combined_custom_concentrations', 'combined_custom_deltaG_error', 'combined_custom_comp_pot']
    
    for param in added_parameters:
        assert param in full_thermo.parameters
    assert 'modelseed_api' in list(full_thermo.variables.keys())           
    
            
    # reinstantiate the model
    bigg_model_path = 'e_coli_core metabolism from BiGG.json'
    model = cobra.io.load_json_model(bigg_model_path)
    
            
def test_build_constraint():
    from scipy.constants import physical_constants, calorie, kilo, R
    Faraday = physical_constants['Faraday constant'][0]#C/mol
    
    # define arbitrary initial conditions to test
    test_reaction = 'PFK'
    model_reaction = full_thermo.model.reactions.get_by_id(test_reaction)
    
    # execute the function    
    built_constraint = full_thermo.build_constraint(object = model_reaction)
      
    # assert results of the function
    for metabolite in full_thermo.model.metabolites:    
        
        # evaluate the thermodynamic calculations
        compartment_potential = 0
        if object.compartment in self.parameters["combined_custom_comp_pot"]:
            compartment_potential = self.parameters["combined_custom_comp_pot"][object.compartment]
            
        msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
        mscpd = ms_api.get_seed_compound(msid)
        constant = mscpd.deltag/calorie + metabolite.charge*Faraday*compartment_potential/kilo
        
        assert type(constant) is float 
        assert constant == built_constraint.ub == built_constraint.lb
    
    import re
    temperature = 298
    for var in built_constraint.variables:
        if re.search('_logconc', str(var)):
            theoretical_coef = -1*R/kilo*temperature
            calculated_coef = float(built_constraint.expression.coeff(var))
            assert theoretical_coef == calculated_coef
    
    # reinstantiate the model
    bigg_model_path = 'e_coli_core metabolism from BiGG.json'
    model = cobra.io.load_json_model(bigg_model_path)
    
    
def test_build_variable():
    
    # define arbitrary initial conditions to test
    test_reaction = 'PFK'
    model_reaction = full_thermo.model.reactions.get_by_id(test_reaction)
    
    # execute the function    
    built_constraint = full_thermo.build_constraint(object = model_reaction)
    
    for metabolite in full_thermo.model.metabolites:   
        full_thermo_conc_var = full_thermo.variables["logconc"][metabolite.id]
        assert full_thermo_conc_var.type == 'continuous'
        assert type(full_thermo_conc_var.lb) is (int or float)
        assert type(full_thermo_conc_var.ub) is (int or float)
        
        full_thermo_dgerr_var = full_thermo.variables["dgerr"][metabolite.id]
        assert full_thermo_dgerr_var.type == 'continuous'
        assert full_thermo_dgerr_var.lb == -full_thermo_dgerr_var.ub