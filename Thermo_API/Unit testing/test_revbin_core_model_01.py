# import the COBRA model
import cobra
bigg_model_path = '..\COBRA function scripts\e_coli_core metabolism from BiGG.json'
model = cobra.io.load_json_model(bigg_model_path)

import optlang

# ------------------------ test the BaseFBA Package ---------------------------------------

# instantiate the class
import sys
package_path = '..\ModelSEEDpy\modelseedpy\\fbapkg'
sys.path.insert(0, package_path)
from revbinpkg import RevBinPkg 

revbin = RevBinPkg(model = model)    
        
def test_init():
    
    # assert results of the model 
    assert type(revbin.model) is cobra.core.model.Model
    assert type(revbin.name) is str
    assert type(revbin.variable_types) is dict
    assert type(revbin.constraint_types) is dict
    assert 'revbin' in list(revbin.variables.keys())
    assert 'revbinR' in list(revbin.constraints.keys())
    assert 'revbinF' in list(revbin.constraints.keys())
    

def test_build_package():
    # build_variable parameters
    lower_bound = 0
    upper_bound = 1
    var_type = "binary"
    constraint_types = {"F": [None, 0],
                        'R': [None, 1000]}
    
    # execute the function
    revbin.build_package()

    # execute the function and assert results of the function
    for reaction in revbin.model.reactions:
        revbin_var = revbin.variables["revbin"][reaction.id]
        assert revbin_var
        assert revbin_var.ub == upper_bound
        assert revbin_var.lb == lower_bound
        assert revbin_var.type == 'binary'
        
        for type in constraint_types:
            revbin_cons = revbin.constraints['revbin{}'.format(type)][reaction.id]
            assert revbin_cons
            assert revbin_cons.lb == constraint_types[type][0]
            assert revbin_cons.ub == constraint_types[type][1]    
            assert revbin_cons.name == '{}_revbin{}'.format(reaction.id, type)