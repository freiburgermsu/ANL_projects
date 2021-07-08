import modelseedpy
from modelseedpy.biochem import from_local
import cobra
import sys

# locally import the ModelSEED API 
import thermodynamic_API_KBase_02
#modelseed_path = '..\\..\\..\\Biofilm growth code\\GSWL code\\ModelSEEDDatabase'
#modelseed = modelseedpy.biochem.from_local(modelseed_path)

# ------------------------ load the test KBase data ---------------------------------------

# general module imports
from optlang.symbolics import Zero, add
from scipy.constants import calorie
from datetime import date
from pandas import DataFrame
import json
import os
import re

# locally import modelseed
import modelseedpy
from modelseedpy.biochem import from_local
modelseed_path = '..\..\..\Biofilm growth code\GSWL code\ModelSEEDDatabase'
modelseed = modelseedpy.biochem.from_local(modelseed_path)

# KBase model import
os.environ["HOME"] = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\cobrakbase'
import cobrakbase
token = '5QQKGJK7BYX7HF7M2TFI3EVJXC7NE67T'
kbase = cobrakbase.KBaseAPI(token)

# remove duplicate compounds from the KBase model 
object_json = kbase.get_object('E_iAH991V2', 93832)
unique_mcs = dict((x['id'], x) for x in object_json['modelcompounds']) 
object_json['modelcompounds'] = list(unique_mcs.values()) # update data without the duplicate id
kbase.save_object('E_iAH991V2', 'freiburgermsu:narrative_1624557251879', 'KBaseFBA.FBAModel', object_json) # saving the object back to KBase with id=E_iAH991V2 to the workspace freiburgermsu:narrative_1624557251879 assigning this type KBaseFBA.FBAModel and with data=object_json

# import the refined KBase individual model
model = kbase.get_from_ws('E_iAH991V2', 93832)

# parse the data
reactions = {}
metabolites = {}
undescribed_compounds = set()
for reaction in object_json['modelreactions']:
    reagents = {}
    reaction_id = reaction['id']
    reaction_name = reaction['name']
    stoich = reaction['modelReactionReagents']
    for reagent in stoich:
        # parse the ModelSEED compound id
        modelseed_id = re.search('(?<=\/)([^\/]+$)', reagent['modelcompound_ref']).group(1)
        
        # parse the ModelSEED compound compartment and refine the ModelSEED ID
        compartment = re.search('(_[a-z][0-9])', modelseed_id).group(1)
        compartment = re.sub('(_)', '', compartment)
        modelseed_id = re.sub('(_[a-z][0-9])', '', modelseed_id)
        
        #print(modelseed_id)
        if re.search('cpd', modelseed_id):
            for key, value in modelseed.get_seed_compound(modelseed_id).data.items():
                if key == 'charge':
                    compound_charge = value
                elif key == 'deltag':
                    compound_gibbs = value / calorie
        
        else:
            print('ERROR: The {} compound is undescribed by the ModelSEED database'.format(modelseed_id))
            undescribed_compounds.add(modelseed_id)
            compound_charge = ''
            compound_gibbs = ''
        
        # create the reagent dictionary
        if modelseed_id not in reagents:
            metabolites[modelseed_id] = {'coefficient': reagent['coefficient'], 
                                      'compartment': compartment,
                                      'charge': compound_charge,
                                      'gibbs (KJ/mol)': compound_gibbs}
            
            reagents[modelseed_id] = metabolites[modelseed_id]

        elif modelseed_id in reagents:
            metabolites[modelseed_id].update({'coefficient_2': reagent['coefficient'], 
                                           'compartment_2': compartment,
                                           'charge_2': compound_charge,
                                           'gibbs_2 (KJ/mol)': compound_gibbs})
            
            reagents[modelseed_id] = metabolites[modelseed_id]
            
    # create the reaction dictionary
    reactions[reaction_id] = {'name': reaction_name,
                              'reagents': reagents}
    
'''with open('2021-06-30_undescribed compounds from the E_iAH991V2 individual model.csv', 'w') as output:
    DataFrame(undescribed_compounds, columns = ['undescribed']).to_csv(output)'''

# review the contents of the data parsing  
'''print(list(reactions.values())[-1])
for reagent, information in reagents.items():
    print(reagent, ': ', information)'''


# ------------------------ test the BaseFBA Package ---------------------------------------

# instantiate the class
package_path = '..\Agronne/ModelSEEDpy/modelseedpy/fbapkg/basefbapkg.py'
sys.path.insert(0, package_path)
import basefbapkg

base = basefbapkg(model = model, name = 'test_model')    

        
def test_init():
    
    # assert results of the model 
    assert base.model is cobrakbase.core.kbasefba.fbamodel.FBAModel
    assert base.name is str
    assert base.variable_types is dict
    assert base.constraint_types is dict
    

def test_validate_parameters():
    # define arbitrary argument content for the function test
    params = {'a':2, 'b':4, 'c':3}
    required = ['a', 'b', 'c']
    defaults = {'a':1, 'b':1, 'c':1, 'd':1}

    # execute the function and assert results of the function
    base.validate_parameters(params, required, defaults)
    assert base.parameters == {'a': 2, 'b': 4, 'c': 3, 'd': 1}
    

def test_build_reaction_variable():
    # define the model instance of the test reaction
    test_reaction = '12ETHDt_c0'
    ub = 133
    lb = 10
    var_type = 'continuous'
    variable_type = 'concentration'
    model_reaction = model.reactions.get_by_id(test_reaction)
    
    # execute the function 
    built_variable = base.build_variable(type = variable_type, lower_bound = lb, upper_bound = ub, vartype = var_type, object = model_reaction)
    
    # assert results of the function
    assert built_variable.ub == ub
    assert built_variable.lb == lb
    assert built_variable.type == var_type
    assert type(built_variable) is tuple
    assert len(base.variables['concentration']) == 1
    assert built_constraint.name == '{}_{}'.format(test_reaction, constraint_type)
    

def test_build_constraint():
    # define arbitrary argument content for the function test
    test_reaction = '12ETHDt_c0'
    ub = 133
    lb = 10
    constraint_type = 'concentration'
    model_reaction = model.reactions.get_by_id(test_reaction)
    
    # execute the function 
    built_constraint = base.build_constraint(type = constraint_type, lower_bound = lb, upper_bound = ub, object = model_reaction)
    
    # assert results of the function
    assert built_constraint.ub == ub
    assert built_constraint.lb == lb
    assert len(base.constraints['concentration']) == 1
    assert built_constraint is optlang.cplex_interface.Constraint
    assert built_constraint.name == '{}_{}'.format(test_reaction, constraint_type)
    
    
def test_all_variables():
    # define the initial conditions of the model
    variables_quantity = 0
    for child in base.childpkgs:
        for kind in child.variables:
            variables_quantity += 1
    for var in base.variables:
        variables_quantity += 1
    
    # execute the function 
    instance_variables = len(base.all_variables())
    
    # assert results of the function
    assert variables_quantity == instance_variables
    
    
def test_all_constraints():
    # define the initial conditions of the model
    constraints_quantity = 0
    for child in base.childpkgs:
        for kind in child.constraints:
            constraints_quantity += 1
    for var in base.constraints:
        constraints_quantity += 1
    
    # execute the function 
    instance_constraints = len(base.all_constraints())
    
    # assert results of the function
    assert constraints_quantity == instance_constraints

    
def test_write_lp_file():   
    # define arbitrary argument content for the function test
    lp_name = 'unit_test_name'

    # execute the function
    base.write_lp_file(lp_name)
    lp_filename = '{}_{}_0.lp'.format(date.today(), export_filename)
    
    # assert results of the function
    assert os.path.exists(lp_name)
    

def clear():
    # define the initial conditions of the model
    variables_quantity = 0
    constraints_quantity = 0
    for type in base.variables:
        for object in base.variables[type]:
            variables_quantity += 1
    for type in base.constraints:
        for object in base.constraints[type]:
            constraints_quantity += 1
            
    # execute the function 
    instance_constraints = base.clear()
    
    # evaluate the reaction names between ModelSEED and COBRA
    assert variables_quantity > len(base.variables)
    assert len(base.variables) == 0
    assert constraints_quantity > len(base.constraints)
    assert len(base.constraints) == 0