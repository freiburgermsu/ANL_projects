# -*- coding: utf-8 -*-

from __future__ import absolute_import
#import logging
import cobra

#Adding a few exception classes to handle different types of errors
class FeasibilityError(Exception):
    """Error in FBA formulation
    """
    if Exception == '':
        pass
    
    if Exception == '':
        pass
    

#Base class for FBA packages
class BaseFBAPkg:
    def __init__(self, model, name, variable_types = {'concentration': 'float', 'lnconc': 'float'}, constraint_types = {'concentration': 'float'}, parent=None):
        '''Intantiate the model
        'model' (COBRA obj): The COBRApy model obj
        'name' (Python obj, string): The name of the model
        'variable_types' (Python obj, dictionary): The types and variables examples for the model variables
        'constraint_types' (Python obj, dictionary): The names and values of the model constraints
        'parent' (Python obj, boolean): The categorical description of the model 
        '''
        import re
        
        self.model = model
        self.name = name
        self.childpkgs = dict()
        self.constraints = dict()
        self.variables = dict()
        self.parameters = dict()
        self.variable_types = variable_types
        self.constraint_types = constraint_types    
        self.thermo_compounds = {}   
        
        BaseFBAPkg.add_variables_and_constraints(self, variable_types = variable_types, constraint_types = constraint_types)

        from scipy.constants import calorie
        import pandas
        import re

        #import modelseedpy      
        #modelseed_path = '..\..\..\Biofilm growth code\GSWL code\ModelSEEDDatabase'
        #modelseed = modelseedpy.biochem.from_local(modelseed_path)
        
        
    def stoichiometry_data(self, obj, view_errors = True, obj_type = 'reaction'):
        ''' Introduce the requisite stoichiometric data to the self instance data
        'obj' (COBRA obj): The name of a COBRA reaction or metabolite, although, the former is the essential intention of the API  
        'obj_type' (Python obj, string): A description of the COBRA obj, which is used to apply the pertinent code for the passed obj
        '''
        import re
        
        # expand the thermodynamic data from the obj argument 
        undescribed_compounds = []
        if obj_type == 'reaction':
            for metabolite in obj.metabolites:
                try:
                    tmfa_name = re.sub('(?i)(_[a-z]$)', '', metabolite.id)
                    self.thermo_compounds[tmfa_name]['stoichiometry'] = obj.metabolites[metabolite]
                    
                except:
                    undescribed_compounds.append(tmfa_name)
                    if view_errors:
                        print('ERROR: The metabolite {} is unrepresented in the thermodynamic database.'. format(tmfa_name))  
        else:
            print('ERROR: The obj_type is not compatible with this API.')
            
        return undescribed_compounds
        
        
    def add_variables_and_constraints(self, variables = {}, variable_types = {}, constraints = {}, constraint_types = {}):
        '''Add arbitrary variables and constraints to a class instance
        'variables' (Python obj, dictionary\list): An iterable set of variables that will be added to the instance veriables
        'constraints' (Python obj, dictionary\list): An iterable set of constraints that will be added to the instance veriables
        'variable_types' (Python obj, dictionary): The types and variables examples for the model variables
        'constraint_types' (Python obj, dictionary): The names and values of the model constraints
        '''
        for variable in variables:
            if variable not in self.variables:
                self.variables[variable] = {}
            
        for kind in variable_types:
            if kind not in self.variables:
                self.variables[kind] = dict()
            if kind not in self.variable_types:
                self.variable_types[kind] = dict()
        
        for constraint in constraints:
            if constraint not in self.constraints:
                self.constraints[constraint] = {}
            
        for kind in constraint_types:
            if kind not in self.constraints:
                self.constraints[kind] = dict()
            if kind not in self.constraint_types:
                self.constraint_types[kind] = dict()
                
                        
    def validate_parameters(self, params, required, defaults):
        '''Validate the passed parameters with the required parameters and update the instance parameters with the default parameters
        'params' (Python obj, dictionary): The parameters and values of the model
        'required' (Python obj, list): The required parameters for the model 
        'defaults' (Python obj, dictionary): The default parameters and values for the model 
        '''
        missing_parameters = set(required) - set(params)
        missing_string = ', '.join(list(missing_parameters))
        if len(missing_parameters) > 1:
            raise ValueError('The required parameters < {} > are missing!'.format(missing_string))
        elif len(missing_parameters) == 1:
            raise ValueError('The required parameter < {} > is missing!'.format(missing_string))
        
        self.parameters = params
        for key in defaults:
            if key not in params:
                self.parameters[key] = defaults[key]
        
        
    def build_variable(self, kind, lower_bound, upper_bound, vartype, obj, obj_type = 'reaction', view_errors = True, errors_quantity = 0):
        '''Create variables of the specified type in the COBRA model
        'kind' (Python obj, string): The variable type within which variables will be created
        'lower_bound' (Python obj, float): The lower bound value for the added variable
        'upper_bound' (Python obj, float): The upper bound value for the added variable
        'vartype' (Python obj, string): The variable type as either 'continuous', 'integer', or 'binary' 
        'obj' (COBRA obj): The COBRA entity into which a variable will be build
        'obj_type' (Python obj): The variable type of the COBRA obj 
        '''
        import re
        
        # assign a variable name based upon the passed arguments
        self.variable_types[kind] = vartype  
        if obj_type == "none":
            count = len(self.variables[kind])
            name = str(count + 1)
        elif obj_type == "string":
            name = obj
        elif obj_type in ['reaction', 'metabolite']:
            name = obj.id
            
        # add an optlang variable, when the variable is undefined
        #raise ValueError('The obj name {} is not recognized by your model'.format(missing_string))
        if name in self.variables[kind]:
            if view_errors:
                print('ERROR: A {} variable already exists, {}.'.format(name, self.variables[kind][name]))
            if not view_errors:
                errors_quantity += 1
            variable_definition = self.variables[kind][name] 
        
        elif name not in self.variables[kind]:
            variable_name = name + "_" + kind
            variable_definition = self.model.problem.Variable(name = variable_name, lb = lower_bound, ub = upper_bound, type = vartype)
            self.model.add_cons_vars(variable_definition)
            self.variables[kind][name] = variable_definition 
            
        return errors_quantity, variable_definition
        
        
    def build_constraint(self, constraint_expression, kind, lower_bound, upper_bound, coef, obj, obj_type = 'reaction', view_errors = True, errors_quantity = 0):
        '''Create constraints for the COBRA model
        'kind' (Python obj, string): The type of the constraint that will be created 
        'lower_bound' (Python obj, float): The lower bound value for the added constraint
        'upper_bound' (Python obj, float): The upper bound value for the added constraint
        'coef' (Python obj, dictionary): The set of coefficients that define the COBRA model
        'obj' (Python obj, string): The variable name when the name is defined
        'obj_type' (Python obj, string): The variable type of the COBRA obj 
        '''
        from optlang.symbolics import Zero
        
        # assign a constraint name based upon the passed arguments
        if obj_type == "none":
            count = len(self.constraints[type])
            name = str(count + 1)
        elif obj_type == "string":
            name = obj
        elif obj_type in ['reaction', 'metabolite']:
            self.constraint_types[kind] = coef
            name = obj.id
                   
        # add an optlang constraint, when the constraint is undefined 
        if name in self.constraints[kind]:
            if view_errors:
                print('ERROR: A {}_{} constraint already exists, {}.'.format(name, kind, self.constraints[kind][name]))
            elif not view_errors:
                errors_quantity += 1
        else:
            constraint_name = '{}_{}'.format(name, kind)
            self.constraints[kind][name] = self.model.problem.Constraint(expression = Zero, lb = lower_bound, ub = upper_bound, name = constraint_name)
            self.model.add_cons_vars(self.constraints[kind][name])
            self.model.solver.update()
            if len(coef) > 0:
                self.constraints[kind][name].set_linear_coefficients(coef)
                
            self.model.solver.update()
        
        return self.constraints[kind][name]
    
    
    def all_variables(self):
        '''Quantify the variables in the class
        '''
        vars = {}
        for child in self.childpkgs:
            for kind in child.variables:
                vars[kind] = child.variables[kind]

        for kind in self.variables:
            vars[kind] = self.variables[kind]
        
        return vars
    
    
    def all_constraints(self):
        '''Quantify the constraints in the class
        '''
        consts = {}
        for child in self.childpkgs:
            for kind in child.constraints:
                consts[kind] = child.constraints[kind]

        for kind in self.constraints:
            consts[kind] = self.constraints[kind]
            
        return consts
    
    
    def revert_to_original(self, cobra_model_path):
        '''Remove changed components of the COBRA model
        'model_path' (Python obj, string): The path string of the COBRA model
        '''               
        global model
        # remove added variables and constants from the model by re-uploading the COBRA model  
        model = cobra.io.load_json_model(cobra_model_path)
        
        return model
    

    def write_lp_file(self, model, export_filename = 'test'):
        '''Export the LP file of the COBRA model
        'model' (COBRA obj): The COBRA model that is expanded through this API
        'export_filename' (Python obj, string): The string of the lp file that will be exported
        '''
        from datetime import date
        import os
        
        import_iteration = 0
        lp_filename = '{}_{}_{}.lp'.format(date.today(), export_filename, import_iteration)
        while os.path.exists(lp_filename):
            import_iteration += 1
            lp_filename = '{}_{}_{}.lp'.format(date.today(), export_filename, import_iteration)
            
        with open(lp_filename, 'w') as lp_file:
            lp_file.write(str(model.solver))
         
        
# ---------------------------------------------- Revbin -------------------------------------------------

# The base class for FBA packages is inherited
class RevBinPkg(BaseFBAPkg):
    def __init__(self, model, obj):
        '''Redefining the inherited __init__ function
        'model' (COBRA obj): The COBRApy FBA model
        '''
        BaseFBAPkg.__init__(self, model = model, name = "reversible binary", variable_types = {"revbin":"reaction", "forv":"reaction", "revv":"reaction"}, constraint_types = {"revbinF":"reaction", "revbinR":"reaction"})
       
        
    def build_constraint(self, obj, view_errors = True, obj_type = 'reaction'):
        '''Build constraints through the inherited function and the calculated coefficient fluxes
        'obj' (Python obj, string): The variable name when the name is defined
        'obj_type' (Python obj, string): The variable type of the COBRA obj 
        '''
        from optlang.symbolics import Zero
        
        # define the constraints of the system
        if obj_type == 'reaction':
            # define the variables that are used in the constraints
            if obj.id not in self.variables['revbin']:
                BaseFBAPkg.build_variable(self, kind = "revbin", lower_bound = 0, upper_bound = 1, vartype = "binary", obj = obj, view_errors = view_errors, obj_type = 'reaction')
                
            revbin_variable = self.variables['revbin'][obj.id] 
            
            coef = {revbin_variable:-1000, obj.forward_variable:1}
            built_forward_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = "revbinF", lower_bound = None, view_errors = view_errors, upper_bound = 0, coef = coef, obj = obj)

            coef = {revbin_variable:1000, obj.reverse_variable:1}
            built_backward_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = "revbinR", lower_bound = None, view_errors = view_errors, upper_bound = 1000, coef = coef, obj = obj)
        
        return built_backward_constraint
    
    
    def build_package(self, variable_kind = 'revbin', filter = None):
        '''Build variables and constraints through the inherited function
        'filter' (Python obj, list): The accepted list of reactions that will be built into the model
        '''
        added_constraints = []
        if filter == None:
            filter = self.variables
            
        for reaction in self.model.reactions:
            # Unfiltered reactions are constructed through the aforementioned functions
            if reaction.id not in filter:
                self.build_constraint(obj = reaction)
                added_constraints.append(reaction.id)
                
        if len(added_constraints) > 0:
            constrained_variables = ', '.join(added_constraints)
            print('The constraints for {} were added to the model.'.format(constrained_variables))
        else:
            print('ERROR: No reactions were added to the model.')
        

                
# ------------------------------------------ Simple Thermo package ------------------------------------------                

# The base class for FBA packages is inherited
class SimpleThermoPkg(RevBinPkg):
    def __init__(self, model, thermo_compounds, view_errors = True):
        '''Redefining the inherited __init__ function
        'model' (COBRA obj): The COBRApy FBA model
        'obj' (COBRA obj): The COBRA reaction\metabolite, or other entity, that should be constrained
        'thermo_compounds' (Python obj, dictionary): The thermodynamic dataset for the model compounds
        '''                               
        # execute the parent __init__ and arbitrarily assign potentials to each metabolite 
        BaseFBAPkg.__init__(self, model = model, name = "simple thermo", variable_types = {"potential":"metabolite", "revbin":"reaction"}, constraint_types = {"simple_thermo":"reaction"})
        
        # inherit the RevBinPkg variables and constraints
        BaseFBAPkg.add_variables_and_constraints(self, variable_types = {"revbin":"reaction", "forv":"reaction", "revv":"reaction"}, constraint_types = {"revbinF":"reaction", "revbinR":"reaction"})
        
        # store the thermodynamic data of the model
        self.thermo_compounds = thermo_compounds
            
            
    def build_constraint(self, obj, obj_type = 'reaction', view_errors = True, errors_quantity = 0, kind = 'simple_thermo'):
        '''Build constraints through the inherited function and the calculated variable coeffiients 
        'obj' (Python obj, string): The variable name when the name is defined
        'obj_type' (Python obj, string): The variable type of the COBRA obj 
        '''
        from optlang.symbolics import Zero
        import re
            
        if obj.id not in self.variables['revbin']:
            BaseFBAPkg.build_variable(self, kind = "revbin", lower_bound = 0, upper_bound = 1, vartype = "binary", obj = obj, view_errors = view_errors, obj_type = 'reaction')
            
        '''different_constraints = set(model.constraints) - set(self.constraints)
        for constraint in different_constraints:
            self.constraints[constraint.name] = constraint.expression'''
        
        BaseFBAPkg.build_variable(self, kind = "potential", lower_bound = 0, upper_bound = 1000, vartype = "continuous", view_errors = view_errors, obj = obj)
        
        constraint_name = '{}_simple_thermo'.format(obj.id)
        coef = {}
        if obj_type == 'reaction':
            if obj.id not in self.constraints['simple_thermo']:
                RevBinPkg.build_constraint(self, obj = obj, obj_type = 'reaction')
                coef[self.variables["revbin"][obj.id]] = 1000
                
                for metabolite in obj.metabolites: 
                    errors_quantity, self.variables["potential"][metabolite.id] = BaseFBAPkg.build_variable(self, kind = "potential", lower_bound = 0, upper_bound = 1000, vartype = "continuous", obj = metabolite, obj_type = 'metabolite', view_errors = view_errors, errors_quantity = errors_quantity)
                    stoichiometry = obj.metabolites[metabolite]
                    coef[self.variables["potential"][metabolite.id]] = stoichiometry
                    
                BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = "simple_thermo", lower_bound = 0, upper_bound = 1000, coef = coef, obj = obj,  view_errors = view_errors, obj_type = 'reaction')
            
            elif obj.id in self.constraints['simple_thermo']:
                if view_errors:
                    print('ERROR: A {} constraint already exists, {}.'.format(constraint_name, self.constraints['simple_thermo'][obj.id]))
                elif not view_errors:
                    errors_quantity += 1      
        
        return errors_quantity

        # calculate the metabolite delta_g in the context of the reaction stoichiometry
        '''delta_g += self.thermo_compounds[tmfa_name]['stoichiometry'] * (self.thermo_compounds[tmfa_name]['gibbs'] + self.variables['concentration_potential'][metabolite.id] + self.variables['electro_potential'][metabolite.id])'''
    
        
        
    def build_package(self, kind = 'potential', view_errors = True, errors_quantity = 0, filter = None):
        '''Build variables and constraints through the inherited function
        'filter' (Python obj, list): The accepted list of reactions that will be built into the model
        '''
        from optlang.symbolics import Zero
        import re
        
        # create thermodynamic constraints and the associated variables
        for reaction in self.model.reactions:   
            errors_quantity = self.build_constraint(obj = reaction, kind = "simple_thermo", obj_type = 'reaction', view_errors = view_errors, errors_quantity = errors_quantity)
            
        if not view_errors:
            print('Quantity of errors: ', errors_quantity)
                
                
# ------------------------------------------ Full Thermo package ------------------------------------------

# The base class for FBA packages is inherited
class FullThermoPkg(SimpleThermoPkg):
    def __init__(self, model, thermo_compounds, thermodynamics_data_type, obj_type = 'reaction'):
        '''Redefining the inherited __init__ function and importing thermodynamic data
        'model' (COBRA obj): The COBRApy FBA model
        'obj' (COBRA obj): The COBRA reaction\metabolite, or other entity, that should be constrained        
        'thermo_reactions' (Python obj, dictionary): The thermodynamic dataset for the model reactions
        'thermo_compounds' (Python obj, dictionary): The thermodynamic dataset for the model compounds
        'thermodynamics_data_type' (Python obj, string): A description of the thermodynamic data which will govern how the data is parsed through the code
        'obj_type' (Python obj, string): The variable type of the COBRA obj 
        '''       
        # execute the base __init__ file
        BaseFBAPkg.__init__(self, model = model, name = "full_thermo", variable_types = {}, constraint_types = {'concentration_potential': 'metabolite', 'electrochemical_potential':'metabolite', 'full_thermo': 'metabolite'})    
        
        # inherit the RevBinPkg variables and constraints
        previous_variables = {"revbin":"reaction", "forv":"reaction", "revv":"reaction", 'ln_concentration': 'metabolite'}
        previous_constraints = {"revbinF":"reaction", "revbinR":"reaction"}
        BaseFBAPkg.add_variables_and_constraints(self, variable_types = previous_variables, constraint_types = previous_constraints)
        
        # parameterize the initial chemical concentrations 
        self.parameters['concentration_potential'] = {}
        self.parameters['electro_potential'] = {}
        self.parameters['total_energy'] = {}
        
        self.validate_parameters(params = self.parameters, required = [], defaults = {
            "default_conc_range": [0.001, 20],  # an arbitrary range
            "custom_concentration_constraints": {"glc-D": [10,10],  #cpd00027,  
                                                'co2': [20, 24], #cpd00011, as bicarbonate, E. B. Brown and Richard L. Clancy, 1964
                                                'h': [0.0053, 0.0053], #cpd00067, Battaglia, Hellegers, & Seeds, 1965
                                                'o2': [0.672, 0.672] #cpd00007, 0.3 mL / dL serum O2 concentration
                                                 },
            'compartment_charge': {'c': 2, # arbitrary value
                                   'e': 0 # by defintion of a zero ph gradient between the metabolite compartment and the extracellular compartment
                                  },
            'activity_coefficient': {"glc-D": 0.94,  # arbitrary value 
                                    'co2': 0.9, # arbitrary value
                                    'h': 0.98, # arbitrary value
                                    'o2': 0.95 # arbitrary value
                                     }
        })
        
        # add the thermodynamic data to this instance
        self.thermo_compounds = thermo_compounds
        self.thermo_data_type = thermodynamics_data_type
        
            
    def build_data(self, obj, coef = {}, view_errors = True, obj_type = 'reaction', errors_quantity = 0):
        ''' The data will be processed into the code instances
        'obj' (COBRA obj): The COBRA reaction\metabolite, or other entity, that should be constrained      
        'thermodynamics_data_type' (Python obj, string): A description of the thermodynamic data which will govern how the data is parsed through the code
        'obj_type' (Python obj, string): The variable type of the COBRA obj 
        '''
        from scipy.constants import physical_constants, kilo, R
        from numpy import log as ln 
        import re
        F = physical_constants['Faraday constant'][0]

        # determine stoichiometric values from the specified dataset
        if self.thermo_data_type == 'dictionary':
            undescribed_compounds = BaseFBAPkg.stoichiometry_data(self, view_errors = view_errors, obj = obj, obj_type = 'reaction')
        
        # calculate the total energy of a reaction based upon the potentials for each metabolite               
        tmfa_name = ''
        if obj_type == 'reaction' and self.thermo_data_type == 'dictionary':
            for metabolite in obj.metabolites:
                tmfa_name = re.sub('(?i)(_[a-z]$)', '', metabolite.id)
                if tmfa_name not in undescribed_compounds:          
                    
                    # determine the concentration parameters for the metabolite
                    if tmfa_name not in self.parameters['custom_concentration_constraints']:
                        activity_coefficient = 1
                        concentrations = self.parameters['default_conc_range']                        
                        #sum(ln(conc) for conc in concentrations) / len(concentration_range)
                    elif tmfa_name in self.parameters['custom_concentration_constraints']:
                        activity_coefficient = self.parameters['activity_coefficient'][tmfa_name]
                        concentrations = self.parameters['custom_concentration_constraints'][tmfa_name]
                        
                    # calculate the concentration potential for the metabolite
                    errors_quantity, self.variables['ln_concentration'][metabolite.id] = BaseFBAPkg.build_variable(self, kind = "ln_concentration", lower_bound = concentrations[0], upper_bound = concentrations[1], vartype = "continuous", obj = metabolite, view_errors = view_errors, errors_quantity = errors_quantity, obj_type = 'metabolite')
                    
                    self.parameters['concentration_potential'][metabolite.id] = self.thermo_compounds[tmfa_name]['stoichiometry'] * self.variables['ln_concentration'][metabolite.id] * ln(activity_coefficient)
                    
                    coef[self.variables['ln_concentration'][metabolite.id]] = self.thermo_compounds[tmfa_name]['stoichiometry'] * ln(activity_coefficient)

                    # calculate the electrochemical potential term
                    psi_electro_potential = 33.33 * self.parameters['compartment_charge'][metabolite.compartment] - 143.33  # millivolts, equation 9 from the TMFA paper 
                    self.parameters['electro_potential'][metabolite.id] = psi_electro_potential * F * self.thermo_compounds[tmfa_name]['charge'] * kilo
            
        return undescribed_compounds
            
    
    def build_constraint(self, obj, filter = 'default_constraints', coef = {}, view_errors = True, errors_quantity = 0, obj_type = 'reaction'): #, modelseed = modelseed):
        '''Build constraints through the inherited function and the calculated variable coeffiients 
        'obj' (Python obj, string): The variable name when the name is defined
        Notes - Equation 14 in the TMFA paper, with the addition of the (charge * compartment_potential) term?
        '''
        from optlang.symbolics import Zero
        from scipy.constants import R
        from numpy import log as ln 
        import re
        
        if obj_type == 'reaction':
            constant = 20 # arbitrary
            if filter == 'default_constraints':
                filter = self.model.constraints
            constraint_name = '{}_full_thermo'.format(obj.id)
            if constraint_name not in filter:
                undescribed_compounds = self.build_data(view_errors = view_errors, coef = coef, obj = obj, obj_type = obj_type, errors_quantity = errors_quantity) 

                # calculate the potential components for the constraint expression calculation   
                temperature = 25 # degrees kelvin  
                if self.thermo_data_type == 'dictionary':
                    for metabolite in obj.metabolites:
                        tmfa_name = re.sub('(?i)(_[a-z]$)', '', metabolite.id)
                        if tmfa_name not in undescribed_compounds:
                            delta_g = self.thermo_compounds[tmfa_name]['gibbs']    
                            electro_potential = self.parameters['electro_potential'][metabolite.id]
                            concentration_potential = self.parameters['concentration_potential'][metabolite.id]            

                            self.parameters['total_energy'][metabolite.id] = delta_g + R * temperature * concentration_potential + electro_potential       

                        else:
                            if view_errors:
                                print('ERROR: The {} compound is not described in the data, thus, the variable was not created.'.format(tmfa_name))

                constraint_names = ['{}_revbinR'.format(obj.id), '{}_revbinF'.format(obj.id)]
                if any(constraint not in self.model.constraints for constraint in constraint_names):
                    RevBinPkg.build_constraint(self, obj = obj, obj_type = 'reaction')
                coef[self.variables["revbin"][obj.id]] = 1000
                built_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = 'full_thermo', lower_bound = constant, upper_bound = constant, coef = coef, obj = obj, view_errors = view_errors, obj_type = 'reaction')
            
            else:
                built_constraint = None
                obj.upper_bound = obj.lower_bound = constant
                print('The {} contraint is updated in the model: {}'.format(constraint_name, self.model.constraints[constraint_name], obj.upper_bound, obj.lower_bound))
                
                if view_errors:
                    print('ERROR: The {}_full_thermo constraint already exists in the model.'.format(obj.id))
            
        else:
            built_constraint = None
            if view_errors:
                print('ERROR: The {} object type is not supported by the API'.format(obj_type))
            
        return built_constraint

        
    def build_package(self, kind = "concentration", view_errors = True, errors_quantity = 0):
        '''Build variables and constraints for a final model package
        'filter' (Python obj, list): The accepted list of reactions that will be built into the model
        Notes - The concentrations are expressed in millimolar
        '''           
        from optlang.symbolics import Zero
        import re

        # create thermodynamic constraints and the associated variables
        for reaction in self.model.reactions: 
            self.build_constraint(obj = reaction, coef = {}, obj_type = 'reaction', view_errors = view_errors, errors_quantity = errors_quantity)
            
        if not view_errors:
            print('Quantity of errors: ', errors_quantity)