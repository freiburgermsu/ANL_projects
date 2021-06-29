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
        'model' (COBRA object): The COBRApy model object
        'name' (Python object, string): The name of the model
        'variable_types' (Python object, dictionary): The types and variables examples for the model variables
        'constraint_types' (Python object, dictionary): The names and values of the model constraints
        'parent' (Python object, boolean): The categorical description of the model 
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
        
        
    def stoichiometry_data(self, object, view_errors = True, object_type = 'reaction'):
        ''' Introduce the requisite stoichiometric data to the self instance data
        'object' (COBRA object): The name of a COBRA reaction or metabolite, although, the former is the essential intention of the API  
        'object_type' (Python object, string): A description of the COBRA object, which is used to apply the pertinent code for the passed object
        '''
        import re
        
        # expand the thermodynamic data from the object argument 
        undescribed_compounds = []
        if object_type == 'reaction':
            for metabolite in object.metabolites:
                try:
                    tmfa_name = re.sub('(?i)(_[a-z]$)', '', metabolite.id)
                    self.thermo_compounds[tmfa_name]['stoichiometry'] = object.metabolites[metabolite]
                    
                except:
                    undescribed_compounds.append(tmfa_name)
                    if view_errors:
                        print('ERROR: The metabolite {} is unrepresented in the thermodynamic database.'. format(tmfa_name))  
        else:
            print('ERROR: The object_type is not compatible with this API.')
            
        return undescribed_compounds
        
        
    def add_variables_and_constraints(self, variables = {}, variable_types = {}, constraints = {}, constraint_types = {}):
        '''Add arbitrary variables and constraints to a class instance
        'variables' (Python object, dictionary\list): An iterable set of variables that will be added to the instance veriables
        'constraints' (Python object, dictionary\list): An iterable set of constraints that will be added to the instance veriables
        'variable_types' (Python object, dictionary): The types and variables examples for the model variables
        'constraint_types' (Python object, dictionary): The names and values of the model constraints
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
        'params' (Python object, dictionary): The parameters and values of the model
        'required' (Python object, list): The required parameters for the model 
        'defaults' (Python object, dictionary): The default parameters and values for the model 
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
        
        
    def build_variable(self, kind, lower_bound, upper_bound, vartype, object, object_type = 'reaction', view_errors = True, errors_quantity = 0):
        '''Create variables of the specified type in the COBRA model
        'kind' (Python object, string): The variable type within which variables will be created
        'lower_bound' (Python object, float): The lower bound value for the added variable
        'upper_bound' (Python object, float): The upper bound value for the added variable
        'vartype' (Python object, string): The variable type as either 'continuous', 'integer', or 'binary' 
        'object' (COBRA object): The COBRA entity into which a variable will be build
        'object_type' (Python object): The variable type of the COBRA object 
        
        '''
        import re
        
        # assign a variable name based upon the passed arguments
        self.variable_types[kind] = vartype  
        if object_type == "none":
            count = len(self.variables[kind])
            name = str(count + 1)
        elif object_type == "string":
            name = object
        elif object_type in ['reaction', 'metabolite']:
            name = object.id
            
        # add an optlang variable, when the variable is undefined
        #raise ValueError('The object name {} is not recognized by your model'.format(missing_string))
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
        
        
    def build_constraint(self, constraint_expression, kind, lower_bound, upper_bound, coef, object, object_type = 'reaction', view_errors = True, errors_quantity = 0):
        '''Create constraints for the COBRA model
        'kind' (Python object, string): The type of the constraint that will be created 
        'lower_bound' (Python object, float): The lower bound value for the added constraint
        'upper_bound' (Python object, float): The upper bound value for the added constraint
        'coef' (Python object, dictionary): The set of coefficients that define the COBRA model
        'object' (Python object, string): The variable name when the name is defined
        'object_type' (Python object, string): The variable type of the COBRA object 
        '''
        from optlang.symbolics import Zero
        
        # assign a constraint name based upon the passed arguments
        if object_type == "none":
            count = len(self.constraints[type])
            name = str(count + 1)
        elif object_type == "string":
            name = object
        elif object_type in ['reaction', 'metabolite']:
            self.constraint_types[kind] = coef
            name = object.id
                   
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
    
    
    def build_data(self, object, thermo_compounds, view_errors = True, object_type = 'reaction', thermodynamics_data_type = 'dictionary'):
        ''' The data will be processed into the code instances
        'object' (COBRA object): The COBRA reaction\metabolite, or other entity, that should be constrained      
        'thermodynamics_data_type' (Python object, string): A description of the thermodynamic data which will govern how the data is parsed through the code
        'object_type' (Python object, string): The variable type of the COBRA object 
        '''
        from scipy.constants import physical_constants, kilo, R
        from numpy import log as ln 
        import re
        F = physical_constants['Faraday constant'][0]

        
        # introduce stoichiometric values for the specified reaction
        self.thermo_compounds = thermo_compounds
        if thermodynamics_data_type == 'dictionary':
            undescribed_compounds = BaseFBAPkg.stoichiometry_data(self, view_errors = view_errors, object = object, object_type = 'reaction')
            
        
        # calculate the total energy of a reaction based upon the values for each constitutent metabolite in the reaction 
        delta_g = 0
        sum_concentration_potential = 0
        sum_electro_potential = 0
        temperature = 25 # degrees kelvin  
        self.variables['concentration_potential'] = {}
        self.variables['electro_potential'] = {}
        self.parameters['activity_coefficient'] = {}
        self.variables['total_energy'] = {}
        tmfa_name = ''
        if object_type == 'reaction' and thermodynamics_data_type == 'dictionary':
            for metabolite in object.metabolites:
                tmfa_name = re.sub('(?i)(_[a-z]$)', '', metabolite.id)
                if tmfa_name not in undescribed_compounds:               
                    # calculate the concentration range for the chemical specie
                    if tmfa_name not in self.parameters['custom_concentration_constraints']:
                        concentration_range = self.parameters['default_conc_range']
                        self.thermo_compounds[tmfa_name]['concentration'] = sum(concentration_range) / len(concentration_range)
                    elif tmfa_name in self.parameters['custom_concentration_constraints']:
                        concentration_range = self.parameters['custom_concentration_constraints'][tmfa_name]
                        self.thermo_compounds[tmfa_name]['concentration'] = sum(concentration_range) / len(concentration_range)            


                    # calculate the electrochemical potential term
                    if metabolite.compartment == 'c':
                        ph_gradient = 2 # arbitrary value
                    elif metabolite.compartment == 'e':
                        ph_gradient = 0 # by defintion of a zero ph gradient between the metabolite compartment and the extracellular compartment

                    psi_electro_potential = 33.33 * ph_gradient - 143.33  # millivolts, equation 9 from the TMFA paper 
                    self.variables['electro_potential'][metabolite.id] = psi_electro_potential * F * self.thermo_compounds[tmfa_name]['charge'] * kilo


                    # calculate the concentration potential term from the average concentration
                    if metabolite.id not in self.parameters['activity_coefficient']:
                        activity_coefficient = 1
                    elif metabolite.id in self.parameters['activity_coefficient']:
                        activity_coefficient = self.parameters['activity_coefficient'][metabolite.id]

                    self.variables['concentration_potential'][metabolite.id] = self.thermo_compounds[tmfa_name]['stoichiometry'] * ln(self.thermo_compounds[tmfa_name]['concentration'] * activity_coefficient)


                    # sum the all energetic descriptions of each metabolite in a reaction
                    sum_concentration_potential += self.variables['concentration_potential'][metabolite.id]
                    sum_electro_potential += self.variables['electro_potential'][metabolite.id]
                    delta_g += self.thermo_compounds[tmfa_name]['gibbs']     
            
            
        # calculation of the total energetic potential of the reaction based upon the metabolite calculations
        self.variables['total_energy'][object.id] = delta_g + R * temperature * sum_concentration_potential + sum_electro_potential  
        
        '''        
        # calculate the total energy of a reaction based upon the values for each constitutent metabolite in the reaction 
        delta_g = 0
        self.parameters['activity_coefficient'] = {}
        for metabolite in object.metabolites:
            tmfa_name = re.sub('(_.)', '', metabolite.id)
            
        if object_type == 'reaction':
            delta_g = 0
            for metabolite in object.metabolites:
                tmfa_name = re.sub('(_.)', '', metabolite.id)
                if thermodynamics_data_type == 'dictionary':
                    delta_g += self.thermo_compounds[tmfa_name]['stoichiometry'] * (self.thermo_compounds[tmfa_name]['gibbs'] + self.constraints['concentration_potential'][metabolite.id] + self.constraints['electro_potential'][metabolite.id])
        '''
        return undescribed_compounds
        
        
        
    '''
    def build_variable(self, object): #, modelseed = modelseed):   
        #flux variability analysis? 
        
        from numpy import log as ln
        if object.id in self.parameters["custom_concentration_constraints"]:
            lb = ln(self.parameters["custom_concentration_constraints"][object.id][0])
            ub = ln(self.parameters["custom_concentration_constraints"][object.id][1])
        else:
            lb = ln(self.parameters["default_min_conc"])
            ub = ln(self.parameters["default_max_conc"])
            
        return BaseFBAPkg.build_variable(self, type = "lnconc", lower_bound = lb, upper_bound = ub, vartype = "continuous", object = object)
    '''
    
    
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
        'model_path' (Python object, string): The path string of the COBRA model
        '''               
        global model
        # remove added variables and constants from the model by re-uploading the COBRA model  
        model = cobra.io.load_json_model(cobra_model_path)
        
        return model
    

    def write_lp_file(self, model, export_filename = 'test'):
        '''Export the LP file of the COBRA model
        'model' (COBRA object): The COBRA model that is expanded through this API
        'export_filename' (Python object, string): The string of the lp file that will be exported
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
    def __init__(self, model, object):
        '''Redefining the inherited __init__ function
        'model' (COBRA object): The COBRApy FBA model
        '''
        BaseFBAPkg.__init__(self, model = model, name = "reversible binary", variable_types = {"revbin":"reaction", "forv":"reaction", "revv":"reaction"}, constraint_types = {"revbinF":"reaction", "revbinR":"reaction"})
       
        
    def build_constraint(self, object, view_errors = True, object_type = 'reaction'):
        '''Build constraints through the inherited function and the calculated coefficient fluxes
        'object' (Python object, string): The variable name when the name is defined
        'object_type' (Python object, string): The variable type of the COBRA object 
        '''
        from optlang.symbolics import Zero
        
        # define the constraints of the system
        if object_type == 'reaction':
            # define the variables that are used in the constraints
            if object.id not in self.variables['revbin']:
                BaseFBAPkg.build_variable(self, kind = "revbin", lower_bound = 0, upper_bound = 1, vartype = "binary", object = object, view_errors = view_errors, object_type = 'reaction')
                
            revbin_variable = self.variables['revbin'][object.id] 
            
            coef = {revbin_variable:-1000, object.forward_variable:1}
            built_forward_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = "revbinF", lower_bound = None, view_errors = view_errors, upper_bound = 0, coef = coef, object = object)

            coef = {revbin_variable:1000, object.reverse_variable:1}
            built_backward_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = "revbinR", lower_bound = None, view_errors = view_errors, upper_bound = 1000, coef = coef, object = object)
        
        return built_backward_constraint
    
    
    def build_package(self, variable_kind = 'revbin', filter = None):
        '''Build variables and constraints through the inherited function
        'filter' (Python object, list): The accepted list of reactions that will be built into the model
        '''
        added_constraints = []
        if filter == None:
            filter = self.variables
            
        for reaction in self.model.reactions:
            # Unfiltered reactions are constructed through the aforementioned functions
            if reaction.id not in filter:
                self.build_constraint(object = reaction)
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
        'model' (COBRA object): The COBRApy FBA model
        'object' (COBRA object): The COBRA reaction\metabolite, or other entity, that should be constrained
        'thermo_compounds' (Python object, dictionary): The thermodynamic dataset for the model compounds
        '''                               
        # execute the parent __init__ and arbitrarily assign potentials to each metabolite 
        BaseFBAPkg.__init__(self, model = model, name = "simple thermo", variable_types = {"potential":"metabolite", "revbin":"reaction"}, constraint_types = {"simple_thermo":"reaction"})
        
        # inherit the RevBinPkg variables and constraints
        BaseFBAPkg.add_variables_and_constraints(self, variable_types = {"revbin":"reaction", "forv":"reaction", "revv":"reaction"}, constraint_types = {"revbinF":"reaction", "revbinR":"reaction"})
        
        # store the thermodynamic data of the model
        self.thermo_compounds = thermo_compounds
            
            
    def build_constraint(self, object, object_type = 'reaction', view_errors = True, errors_quantity = 0, kind = "potential"):
        '''Build constraints through the inherited function and the calculated variable coeffiients 
        'object' (Python object, string): The variable name when the name is defined
        'object_type' (Python object, string): The variable type of the COBRA object 
        '''
        from optlang.symbolics import Zero
        import re
            
        if object.id not in self.variables['revbin']:
            BaseFBAPkg.build_variable(self, kind = "revbin", lower_bound = 0, upper_bound = 1, vartype = "binary", object = object, view_errors = view_errors, object_type = 'reaction')
            
        if object_type == 'reaction':
            if object.reversibility:
                binary = 0
            elif not object.reversibility:
                binary = 1
            else:
                print('ERROR: The reaction object possesses unpredictable data structure.')
            
        '''different_constraints = set(model.constraints) - set(self.constraints)
        for constraint in different_constraints:
            self.constraints[constraint.name] = constraint.expression'''
        
        BaseFBAPkg.build_variable(self, kind = kind, lower_bound = 0, upper_bound = 1000, vartype = "continuous", view_errors = view_errors, object = object)
        
        constraint_name = '{}_simple_thermo'.format(object.id)
        coef = {}
        if object_type == 'reaction':
            if object.id not in self.constraints['simple_thermo']:
                RevBinPkg.build_constraint(self, object = object, object_type = 'reaction')
                coef[self.variables["revbin"][object.id]] = 1000
                
                for metabolite in object.metabolites: 
                    errors_quantity, self.variables[kind][metabolite.id] = BaseFBAPkg.build_variable(self, kind = kind, lower_bound = 0, upper_bound = 1000, vartype = "continuous",  view_errors = view_errors, object = metabolite, object_type = 'metabolite')
                    stoichiometry = object.metabolites[metabolite]
                    coef[self.variables[kind][metabolite.id]] = stoichiometry
                    
                BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = "simple_thermo", lower_bound = 0, upper_bound = 1000, coef = coef, object = object,  view_errors = view_errors, object_type = 'reaction')
            
            elif object.id in self.constraints['simple_thermo']:
                if view_errors:
                    print('ERROR: A {} constraint already exists, {}.'.format(constraint_name, self.constraints['simple_thermo'][object.id]))
                elif not view_errors:
                    errors_quantity += 1
        
        print(errors_quantity)            
        
        
    def build_package(self, kind = 'potential', view_errors = True, errors_quantity = 0, filter = None):
        '''Build variables and constraints through the inherited function
        'filter' (Python object, list): The accepted list of reactions that will be built into the model
        '''
        from optlang.symbolics import Zero
        import re
        
        #RevBinPkg.build_package(self, variable_kind = 'simple_thermo', filter = filter)
        
        # create thermodynamic constraints and the associated variables
        for reaction in self.model.reactions: 
            coef = {}
            for metabolite in reaction.metabolites:
                #tmfa_name = re.sub('(_.)', '', metabolite.id)
                
                errors_quantity, self.variables[kind][metabolite.id] = BaseFBAPkg.build_variable(self, kind = kind, lower_bound = 0, upper_bound = 1000, vartype = "continuous", object = metabolite, object_type = 'metabolite', view_errors = view_errors, errors_quantity = errors_quantity)
                stoichiometry = reaction.metabolites[metabolite]
                coef[self.variables[kind][metabolite.id]] = stoichiometry
                
            BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = "simple_thermo", lower_bound = 0, upper_bound = 1000, coef = coef, object = reaction, view_errors = view_errors, object_type = 'reaction')
            
        if not view_errors:
            print('Quantity of errors: ', errors_quantity)
                
                
# ------------------------------------------ Full Thermo package ------------------------------------------

# The base class for FBA packages is inherited
class FullThermoPkg(SimpleThermoPkg):
    def __init__(self, model, thermodynamics_data_type, object_type = 'reaction'):
        '''Redefining the inherited __init__ function and importing thermodynamic data
        'model' (COBRA object): The COBRApy FBA model
        'object' (COBRA object): The COBRA reaction\metabolite, or other entity, that should be constrained        
        'thermo_reactions' (Python object, dictionary): The thermodynamic dataset for the model reactions
        'thermo_compounds' (Python object, dictionary): The thermodynamic dataset for the model compounds
        'thermodynamics_data_type' (Python object, string): A description of the thermodynamic data which will govern how the data is parsed through the code
        'object_type' (Python object, string): The variable type of the COBRA object 
        '''       
        # execute the base __init__ file
        BaseFBAPkg.__init__(self, model = model, name = "full_thermo", variable_types = {"lnconc":"metabolite", 'full_thermo': 'reactions', "revbin":"reaction"}, constraint_types = {'concentration_potential': 'metabolite', 'electrochemical_potential':'metabolite', 'thermo': 'metabolite', "full_thermo":"reaction"})    
        
        # inherit the RevBinPkg variables and constraints
        BaseFBAPkg.add_variables_and_constraints(self, variable_types = {"revbin":"reaction", "forv":"reaction", "revv":"reaction"}, constraint_types = {"revbinF":"reaction", "revbinR":"reaction"})
        
        # parameterize the initial chemical concentrations 
        self.validate_parameters(params = self.parameters, required = [], defaults = {
            "default_conc_range": [0.001, 20],  # an arbitrary range
            "custom_concentration_constraints": {"glc-D": [10,10],  #cpd00027,  
                                                'co2': [20, 24], #cpd00011, as bicarbonate, E. B. Brown and Richard L. Clancy, 1964
                                                'h': [0.0053, 0.0053], #cpd00067, Battaglia, Hellegers, & Seeds, 1965
                                                'o2': [0.672, 0.672] #cpd00007, 0.3 mL / dL serum O2 concentration
                                                 }
        })
               
    
    def build_constraint(self, object, thermo_compounds, thermodynamics_data_type, view_errors = True, kind = 'full_thermo', object_type = 'reaction'): #, modelseed = modelseed):
        '''Build constraints through the inherited function and the calculated variable coeffiients 
        'object' (Python object, string): The variable name when the name is defined
        Notes - Equation 14 in the TMFA paper, with the addition of the (charge * compartment_potential) term?
        '''
        from optlang.symbolics import Zero
        import re
        
        undescribed_compounds = BaseFBAPkg.build_data(self, view_errors = view_errors, object = object, object_type = object_type, thermodynamics_data_type = thermodynamics_data_type, thermo_compounds = thermo_compounds) 
        
        # calculate the parameters for the constraint expression calculation               
        constant = 20  # arbitrary value
        
        delta_g = 0
        if thermodynamics_data_type == 'dictionary':
            for metabolite in object.metabolites:
                tmfa_name = re.sub('(?i)(_[a-z]$)', '', metabolite.id)
                if tmfa_name not in undescribed_compounds: 
                    delta_g += self.thermo_compounds[tmfa_name]['stoichiometry'] * (self.thermo_compounds[tmfa_name]['gibbs'] + self.variables['concentration_potential'][metabolite.id] + self.variables['electro_potential'][metabolite.id])
            
        if object.reversibility:
            binary = 0
        elif not object.reversibility:
            binary = 1
        else:
            if view_errors:
                print('ERROR: The reaction object possesses unpredictable data structure.')
        
        constraint_names = []
        for constraint in self.model.constraints:
            constraint_names.append(constraint.name)
        
        constraint_name = '{}_{}'.format(object.id, kind)
        coef = {}
        if constraint_name not in constraint_names:
            RevBinPkg.build_constraint(self, object = object, object_type = 'reaction')
            coef[self.variables["revbin"][object.id]] = 1000
            built_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = kind, lower_bound = constant, upper_bound = constant, coef = coef, object = object, view_errors = view_errors, object_type = 'reaction')
            
        else:
            if object_type == 'reaction':
                object.upper_bound = constant
                object.lower_bound = constant
                #model.constraints[constraint_name] = constraint_expression
                if view_errors:
                    print('The {} contraint is updated in the model: {}, ub: {}, lb: {}'.format(constraint_name, model.constraints[constraint_name], object.upper_bound, object.lower_bound))

            built_constraint = None
            
        return built_constraint

        
    def build_package(self, thermo_compounds, kind = "full_thermo", view_errors = True, errors_quantity = 0):
        '''Build variables and constraints for a final model package
        'filter' (Python object, list): The accepted list of reactions that will be built into the model
        Notes - The concentrations are expressed in millimolar
        '''           
        from optlang.symbolics import Zero
        import re

        filter = self.constraints

        # create thermodynamic constraints and the associated variables
        for reaction in self.model.reactions: 
            coef = {}
            for metabolite in reaction.metabolites:
                #tmfa_name = re.sub('(_.)', '', metabolite.id)
                
                errors_quantity, self.variables[kind][metabolite.id] = BaseFBAPkg.build_variable(self, kind = kind, lower_bound = 0, upper_bound = 1000, vartype = "continuous", object = metabolite, object_type = 'metabolite', view_errors = view_errors, errors_quantity = errors_quantity)
                stoichiometry = reaction.metabolites[metabolite]
                coef[self.variables[kind][metabolite.id]] = stoichiometry
                
            BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = kind, lower_bound = 0, upper_bound = 1000, coef = coef, view_errors = view_errors, object = reaction, object_type = 'reaction')
            
        
        # The concentration variable and potential constraint are built
        for metabolite in self.model.metabolites:
            BaseFBAPkg.build_variable(self, kind = kind, lower_bound = 0, upper_bound = 1000, vartype = "continuous", view_errors = view_errors, object = metabolite)
            
        for reaction in self.model.reactions:
            # Unfiltered reactions are constructed through the aforementioned functions
            if reaction.id not in filter:
                self.build_constraint(object = reaction, thermo_compounds = thermo_compounds, view_errors = view_errors, thermodynamics_data_type = 'dictionary')
                
        if not view_errors:
            print('Quantity of errors: ', errors_quantity)