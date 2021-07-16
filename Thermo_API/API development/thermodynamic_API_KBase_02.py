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
        'model' (COBRA or kbase obj): The genome-scale model obj
        'name' (Python obj, string): The name of the model
        'variable_types' (Python obj, dictionary): The types and variables examples for the model variables
        'constraint_types' (Python obj, dictionary): The names and values of the model constraints
        'parent' (Python obj, boolean): The categorical description of the model 
        '''        
        self.model = model
        self.name = name
        self.childpkgs = dict()
        self.constraints = dict()
        self.variables = dict()
        self.parameters = dict()
        self.variable_types = variable_types
        self.constraint_types = constraint_types 
        
        BaseFBAPkg.add_variables_and_constraints(self, variable_types = variable_types, constraint_types = constraint_types)
        
        
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
            print(coef)
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
    
    
    def revert_to_original(self, cobra_model_path = None, kbase = None, kbase_model_id = 'E_iAH991V2', workspace_id = 93832):
        '''Remove alterations of the genomen-scale model
        'model_path' (Python obj, string): The path string of the COBRA model
        Note - either the < cobra_model_path > argument or the < kbase_model_id > and the < workspace_id > arguments must be passed
        '''               
        global model
        
        if cobra_model_path != None:
            model = cobra.io.load_json_model(cobra_model_path)
        elif cobra_model_path == None:
            model = kbase.get_from_ws(kbase_model_id, workspace_id)
            
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
    
    
    def build_package(self, view_errors = True, variable_kind = 'revbin', filter = None):
        '''Build variables and constraints through the inherited function
        'filter' (Python obj, list): The accepted list of reactions that will be built into the model
        '''
        added_constraints = []
        if filter == None:
            filter = self.constraints
            
        for reaction in self.model.reactions:
            # Unfiltered reactions are constructed through the aforementioned functions
            if reaction.id in filter:
                if view_errors:
                    constraint_name = '{}_revbinF'.format(reaction.id)
                    print('ERROR: A {} constraint already exists:\n\t{}\n{}.'.format(constraint_name, self.constraints['revbinF'][reaction.id], self.constraints['revbinR'][reaction.id]))
                '''self.model.solver.remove(reaction.id)'''
                
            self.build_constraint(obj = reaction)
            added_constraints.append(reaction.id)
                
        if len(added_constraints) > 0:
            constrained_variables = ', '.join(added_constraints)
            print('The constraints for {} were added to the model.'.format(constrained_variables))
        else:
            print('ERROR: No reactions were added to the model.')
        

                
# ------------------------------------------ Simple Thermo package ------------------------------------------                
class SimpleThermoPkg(RevBinPkg):
    def __init__(self, model, thermo_compounds, thermodynamics_data_type, view_errors = True):
        '''Redefining the inherited __init__ function
        'model' (COBRA obj): The COBRApy FBA model
        'obj' (COBRA obj): The COBRA reaction\metabolite, or other entity, that should be constrained
        'thermo_compounds' (Python obj, dictionary): The thermodynamic dataset for the model compounds
        '''                               
        # execute the parent __init__ and arbitrarily assign potentials to each metabolite 
        BaseFBAPkg.__init__(self, model = model, name = "simple thermo", variable_types = {"potential":"metabolite", "revbin":"reaction"}, constraint_types = {"simple_thermo":"reaction"})
        
        # inherit the RevBinPkg variables and constraints
        BaseFBAPkg.add_variables_and_constraints(self, variable_types = {"revbin":"reaction", "forv":"reaction", "revv":"reaction"}, constraint_types = {"revbinF":"reaction", "revbinR":"reaction"})
            
        # add the thermodynamic data to this instance
        self.thermo_compounds = thermo_compounds
        self.thermo_data_type = thermodynamics_data_type
        
            
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
                    # the metabolite potential variable is created
                    errors_quantity, self.variables["potential"][metabolite.id] = BaseFBAPkg.build_variable(self, kind = "potential", lower_bound = 0, upper_bound = 1000, vartype = "continuous", obj = metabolite, obj_type = 'metabolite', view_errors = view_errors, errors_quantity = errors_quantity)

                    # the metabolite gibbs potential is defined 
                    if self.thermo_data_type == 'dictionary':
                        tmfa_name = re.sub('(?i)(_[a-z]$)', '', metabolite.id) 
                        try:
                            delta_g = self.thermo_compounds[tmfa_name]['gibbs']   
                        except:
                            delta_g = 0
                            if view_errors:
                                print('ERROR: The {} metabolite is undescribed by the {} dataset'.format(metabolite.id, self.thermo_data_type))
                                
                    elif self.thermo_data_type == 'kbase':
                        try:
                            delta_g = self.thermo_compounds[metabolite.id]['gibbs (KJ/mol)']   
                        except:
                            delta_g = 0
                            if view_errors:
                                print('ERROR: The {} metabolite is undescribed by the {} dataset'.format(metabolite.id, self.thermo_data_type))
                                
                    # the potential variable coefficient is defined
                    stoichiometry = obj.metabolites[metabolite]
                    coef[self.variables["potential"][metabolite.id]] = stoichiometry * delta_g
                    
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
            errors_quantity = SimpleThermoPkg.build_constraint(self, obj = reaction, kind = "simple_thermo", obj_type = 'reaction', view_errors = view_errors, errors_quantity = errors_quantity)
            
        if not view_errors:
            print('Quantity of errors: ', errors_quantity)
                
                
# ------------------------------------------ Full Thermo package ------------------------------------------
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
        previous_variables = {'ln_concentration': 'metabolite', 'potential':'metabolite'}
        previous_constraints = {'simple_thermo':'reaction'}
        BaseFBAPkg.add_variables_and_constraints(self, variable_types = previous_variables, constraint_types = previous_constraints)
        
        # parameterize the initial chemical concentrations 
        self.parameters['concentration_potential'] = {}
        self.parameters['electro_potential'] = {}
        self.parameters['total_energy'] = {}
        
        self.validate_parameters(params = self.parameters, required = [], defaults = {
            "default_conc_range": [0.001, 20],  # an arbitrary range
            "custom_concentration_constraints": {"glc__D": [10,10],  #cpd00027,  
                                                'co2': [20, 24], #cpd00011, as bicarbonate, E. B. Brown and Richard L. Clancy, 1964
                                                'h': [0.0053, 0.0053], #cpd00067, Battaglia, Hellegers, & Seeds, 1965
                                                'o2': [0.672, 0.672] #cpd00007, 0.3 mL / dL serum O2 concentration
                                                 },
            'compartment_charge': {'c': 2, # arbitrary value
                                   'c0': 2, # arbitrary value
                                   'e': 0, # by defintion of a zero ph gradient between the metabolite compartment and the extracellular compartment
                                   'e0': 0 # by defintion of a zero ph gradient between the metabolite compartment and the extracellular compartment
                                  },
            'activity_coefficient': {"glc__D": 0.94,  # arbitrary value 
                                    'co2': 0.9, # arbitrary value
                                    'h': 0.98, # arbitrary value
                                    'o2': 0.95 # arbitrary value
                                     }
        })
        
        # add the thermodynamic data to this instance
        self.thermo_compounds = thermo_compounds
        self.thermo_data_type = thermodynamics_data_type
        
        
    def build_constraint(self, metabolite_obj, psi_electro_compartment_potential, filter = 'default_constraints', coef = {}, view_errors = True, errors_quantity = 0): #, modelseed = modelseed):
        '''Build constraints through the inherited function and the calculated variable coeffiients 
        'metabolite_obj' (COBRA metabolite obj): The COBRA model metabolite upon which constraints are added
        'psi_electro_compartment_potential' (Python obj, float): The electrochemical potential of the specified metabolite in the respective compartment is defined by the user.
        Notes - Equation 14 in the TMFA paper, with the addition of the (charge * compartment_potential) term?
        '''
        from optlang.symbolics import Zero
        from scipy.constants import physical_constants, kilo, R
        from numpy import log as ln 
        import re
        
        # the initial parameters are established 
        F = physical_constants['Faraday constant'][0]
        temperature = 25 # degrees kelvin  
        if filter == 'default_constraints':
            filter = self.model.constraints
            
        # remove an existing constraint 
        constraint_name = '{}_full_thermo'.format(metabolite_obj.id)
        if constraint_name in filter:           
            if view_errors:
                print('ERROR: The {}_full_thermo constraint already exists in the model.'.format(metabolite_obj.id))
            '''self.solver.remove(constraint_name)'''
        
        # calculate the total energy of a reaction based upon the potentials for each metabolite               
        tmfa_name = re.sub('(?i)(_[a-z]$)', '', metabolite_obj.id)

        # determine the metabolite concentration parameters 
        if tmfa_name not in self.parameters['custom_concentration_constraints']:
            activity_coefficient = 1
            concentrations = self.parameters['default_conc_range']                        
            #sum(ln(conc) for conc in concentrations) / len(concentration_range)
        elif tmfa_name in self.parameters['custom_concentration_constraints']:
            activity_coefficient = ln(self.parameters['activity_coefficient'][tmfa_name])
            concentrations = self.parameters['custom_concentration_constraints'][tmfa_name]

        # calculate the metabolite concentration potential 
        errors_quantity, self.variables['ln_concentration'][metabolite_obj.id] = BaseFBAPkg.build_variable(self, kind = "ln_concentration", lower_bound = concentrations[0], upper_bound = concentrations[1], vartype = "continuous", obj = metabolite_obj, view_errors = view_errors, errors_quantity = errors_quantity, obj_type = 'metabolite')
        coef[self.variables['ln_concentration'][metabolite_obj.id]] = activity_coefficient * R * temperature

        # calculate the electrochemical potential term
        #psi_electro_potential = 33.33 * self.parameters['compartment_charge'][metabolite_obj.compartment] - 143.33  ;  in millivolts, equation 9 from the TMFA paper 
        self.parameters['electro_potential'][metabolite_obj.id] = psi_electro_compartment_potential * F * metabolite_obj.charge * kilo

        # calculate the potential components for the constraint expression calculation   
        if self.thermo_data_type == 'dictionary':
            tmfa_name = re.sub('(?i)(_[a-z]$)', '', metabolite_obj.id) 
            try:
                delta_g = self.thermo_compounds[tmfa_name]['gibbs']   
            except:
                delta_g = 0
                if view_errors:
                    print('ERROR: The {} metabolite is undescribed by the {} dataset'.format(metabolite_obj.id, self.thermo_data_type))
                    
        elif self.thermo_data_type == 'kbase':
            try:
                delta_g = self.thermo_compounds[metabolite_obj.id]['gibbs (KJ/mol)']   
            except:
                delta_g = 0
                if view_errors:
                    print('ERROR: The {} metabolite is undescribed by the {} dataset'.format(metabolite_obj.id, self.thermo_data_type))
                    
        # create the constraint
        constant = delta_g + self.parameters['electro_potential'][metabolite_obj.id]
        built_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = 'full_thermo', lower_bound = constant, upper_bound = constant, coef = coef, obj = metabolite_obj, view_errors = view_errors, obj_type = 'metabolite')

        return built_constraint

        
    def build_package(self, electro_compartment_potential_dict, kind = "concentration", view_errors = True, errors_quantity = 0):
        '''Build variables and constraints for a final model package
        'filter' (Python obj, list): The accepted list of reactions that will be built into the model
        Notes - The concentrations are expressed in millimolar
        '''           
        from optlang.symbolics import Zero
        import re

        # create thermodynamic constraints and the associated variables
        for metabolite in self.model.metabolites: 
            if self.thermo_data_type == 'dictionary':
                compartment = metabolite.compartment
                compartment_potential = electro_compartment_potential_dict[compartment]
            elif self.thermo_data_type == 'kbase':
                compartment = re.sub('([0-9])', '', metabolite.compartment)
                compartment_potential = electro_compartment_potential_dict[compartment]
            self.build_constraint(metabolite_obj = metabolite, coef = {}, psi_electro_compartment_potential = compartment_potential, view_errors = view_errors, errors_quantity = errors_quantity)
            
        if not view_errors:
            print('Quantity of errors: ', errors_quantity)