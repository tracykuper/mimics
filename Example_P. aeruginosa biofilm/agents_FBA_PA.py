## IMPORT PYTHON PACKAGES
from multiprocessing import Pool
import sys
import os
import numpy as np
import cobra
import pandas as pd
import math
from py4j.java_gateway import JavaGateway # PY4J package
import time

## INITIALIZE RESULT ARRAY FOR EACH AGENT
full_result = {}

## FUNCTION TO INITIALIZE MULTIPROCESSING WORKERS
def init():
    try:
        print('Multiprocessing worker initialized', flush=True)
    except Exception: # RAISE EXCEPTION
        sarray = gateway.new_array(gateway.jvm.java.lang.String,2) ## STRING FOR PRINTING FROM JAVA
        sarray[0]="MULTIPROCESSING WORKER FAILED"
        gateway.entry_point.Print_Phrase(sarray[0]) ## PRINT EXCEPTION
        sys.stdout.flush()

## FUNCTION TO CONSTAIN AND OPTIMIZE AN AGENT'S METABOLIC MODEL
def run_GENRE(old_biomass,max_biomass,metabolite_concentrations,index,models, metabolic_state,rxns,metabolite_ids,v_patch,dt_rxn,dt_growth):

    model = models[metabolic_state] ## ASSIGN METABOLIC MODEL STATE TO AGENT

    ## SET OUTPUT VALUES FOR AGENT (IN CASE EXCEPTION IS RAISED)
    updated_index = index
    new_biomass = old_biomass
    updated_metabolite_concentrations = metabolite_concentrations
    growth_rate = 0.0
    rxn_flux_all = np.zeros(len(rxns))
    error = ''

    try:
        if old_biomass > 0 and old_biomass<max_biomass: # CHECK AGENT BIOMASS IS ABOVE ZERO
            conversion_metabolite = v_patch * 1/dt_rxn * 1/old_biomass # METABOLITE CONVERSION TERM BETWEEN ABM (mM) <--> GENRE (mmol/ gdWt hr)
            ## SET BOUNDS FOR EXCHANGE METABOLITE REACTIONS
            for ex_rxn,m in zip(metabolite_ids,range(len(metabolite_ids))):
                if ex_rxn in model.exchanges: ## CHECK EXCHANGE METABOLITE REACTION IS IN METABOLIC MODEL
                    if ex_rxn != 'EX_cpd00418_e': ## FOR METABOLITES THAT ARE NOT NITRIC OXIDE
                        model.reactions.get_by_id(ex_rxn).upper_bound = 1000
                        model.reactions.get_by_id(ex_rxn).lower_bound = -metabolite_concentrations[m] * conversion_metabolite;
                    if (ex_rxn == 'EX_cpd00418_e') & (metabolic_state !=2): ## NITRIC OXIDE: CONSTRAIN UPPER BOUND TO SIMULATE FORCED UPTAKE INTO CELL
                        model.reactions.get_by_id(ex_rxn).lower_bound = -1000;
                        model.reactions.get_by_id(ex_rxn).upper_bound = -metabolite_concentrations[m] * conversion_metabolite;
                    if (ex_rxn == 'EX_cpd00418_e') & (metabolic_state ==2): ## NITRIC OXIDE: FOR NO PRODUCERS
                        model.reactions.get_by_id(ex_rxn).upper_bound = 1000
                        model.reactions.get_by_id(ex_rxn).lower_bound = -1000;

            sol = model.optimize() # OPTIMIZE METABOLIC MODEL

            if  sol.objective_value > 1e-2: ## CHECK METABOLIC MODEL SOLUTION IS OPTIMAL

                ## CALCULATE NEW BIOMASS VALUE FOR THIS AGENT
                growth_rate = sol.objective_value # OBTAIN AGENT GROWTH RATE (UNITS: 1/hr)
                new_biomass = old_biomass*math.exp(growth_rate*dt_growth/60) # CALCULATE UPDATED AGENT BIOMASS (UNITS: GRAMS)

                ## CALCULATE NEW METABOLITE CONCENTRATIONS AT THE AGENT'S LOCATION
                for ex_rxn,m in zip(metabolite_ids,range(len(metabolite_ids))):
                    if ex_rxn in model.exchanges:
                        flux = sol.fluxes.loc[ex_rxn]
                        if abs(flux) > 1e-12:
                            updated_metabolite_concentrations[m] = metabolite_concentrations[m]+flux/conversion_metabolite
                        else:
                            updated_metabolite_concentrations[m] = metabolite_concentrations[m]
                        if updated_metabolite_concentrations[m] < 0:
                            updated_metabolite_concentrations[m] = 0.0

            ## SAVE REACTION FLUX VALUES FOR THIS AGENT
            rxn_flux_all = []
            for r in rxns: 
                try: 
                    rxn_flux = model.reactions.get_by_id(r).flux
                    rxn_flux_all.append(rxn_flux)
                except Exception: # RAISE EXCEPTION
                    rxn_flux = 0.0
                    rxn_flux_all.append(rxn_flux)

    except Exception as e: ## IF GENRE OPTIMIZATION FUNCTION FAILS RAISE EXCEPTION
        updated_index = index
        new_biomass = old_biomass
        updated_metabolite_concentrations = metabolite_concentrations
        growth_rate = 0.0
        error = str(e)

    return new_biomass, tuple(updated_metabolite_concentrations), updated_index,growth_rate, tuple(rxn_flux_all),error
    
## FUNCTION TO ACCUMULATE MULTIPROCESSING RUN_GENRE RESULTS
def accumulateResults(result):
    try:
        full_result[result[0]] = result[1]
    except Exception as e:
        sarray = gateway.new_array(gateway.jvm.java.lang.String,2)
        sarray[0]="ACCUMULATE RESULTS FAILED"
        gateway.entry_point.Print_Phrase(sarray[0])
        sarray[0]=str(e)
        gateway.entry_point.Print_Phrase(sarray[0])
        sys.stdout.flush()

## FUNCTION TO RUN MIMICS FOR MULTIPLE SIMULATION TIME STEPS
def run_MIMICS(ncpus, num_dt,models,media,rxns, job_num,num_met_gas,initial_gas_concentrations, D_gas, num_gas_step, num_met_carbon,initial_carbon_concentrations, D_carbon,num_carbon_step, metabolite_ids, v_patch, dt_rxn, dt_growth,initial_biomass,max_biomass,dead_state,lag_phase,output_dir):

    pool = Pool(initializer = init, processes=ncpus) ## INITIALIZE MULTIPROCESSING POOLS

    ## CALCULATE METABOLITE CONVERSION TERM BETWEEN ABM (mM) <--> GENRE (mmol/ gdWt hr)
    dt_rxn = dt_rxn/60; # MINUTES
    dt_rxn = dt_rxn/60; # HOURS
    half_biomass = (max_biomass+initial_biomass)/2 ## MID-POINT BIOMASS VALUE
    conversion_metabolite = v_patch * 1/dt_rxn * 1/half_biomass # METABOLITE CONVERSION TERM BETWEEN ABM (mM) <--> GENRE (mmol/ gdWt hr)

   # CONSTRAIN METABOLIC MODELS ON NUTRIENT MEDIA
    for model_m in models: ## ITERATE OF METABOLIC MODEL STATES
        medium = {} ## DEFINE METABOLITE MEDIUM
        for ex_rxn,i in zip(media['Metabolite ID'],media.index):  ## ITERATE OVER EXCHANGE METABOLITES
            if ex_rxn in model_m.exchanges: ## CHECK EXCHANGE METABOLITE IN METABOLIC MODEL
                food = media.loc[i]['Metabolite concentration (mM)']
                model_m.reactions.get_by_id(ex_rxn).upper_bound = 1000
                model_m.reactions.get_by_id(ex_rxn).lower_bound = -food*conversion_metabolite
                medium[ex_rxn] = food * conversion_metabolite;
        model_m.medium = medium

    ## RUN MIMICS FOR DESIRED NUMBER OF SIMULATION TIME STEPS (num_dt)
    for t in range(num_dt):
        gateway.entry_point.run_ABM(initial_biomass,max_biomass,int(dead_state),lag_phase,int(t)) ## RUN AGENT-BASED MODEL

        ## PRINT SIMULATION TIME STEP FROM JAVA
        sarray = gateway.new_array(gateway.jvm.java.lang.String,2) ## STRING JAVA ARRAY USE TO OUTPUT PRINT STATEMENTS FROM PY4J SERVER
        sarray[0]="Time Step = " + str(t)
        gateway.entry_point.Print_Phrase(sarray[0])

        ## GET VALUES OF AGENT ATTRIBUTES FROM JAVA
        biomass_values = list(gateway.entry_point.getBiomassFromHal()) ## GET BIOMASS OF EACH FROM ABM
        biomass_values = [x / 1e14 for x in biomass_values] ## CONVERT BIOMASS TO GRAMS
        index_values = list(gateway.entry_point.getIndexFromHal()) ## GET INDEX OF EACH AGENT FROM ABM
        patch_values_from_java = np.array(gateway.entry_point.getPatchFromHal_All(int(num_met_gas),int(num_met_carbon))) ## GET METABOLITE CONCENTRATIONS AT EACH AGENT'S LOCATION
        metabolic_states = list(gateway.entry_point.getMetabolicStateFromHal()) ## GET METABOLIC STATE OF EACH AGENT FROM ABM

        ## FORMAT METABOLITE CONCENTRATIONS
        metabolite_concentrations = patch_values_from_java.tolist()

        ## GET ATTRIBUTES OF ALIVE CELLS
        alive_cell_indices = [index for index, value in enumerate(metabolic_states) if value != dead_state]
        index_values = [index_values[index] for index in alive_cell_indices]
        biomass_values = [biomass_values[index] for index in alive_cell_indices]
        metabolic_states = [metabolic_states[index] for index in alive_cell_indices]
        metabolite_concentrations = [metabolite_concentrations[index] for index in alive_cell_indices]
        pos_count = len(biomass_values)  ## GET NUMBER OF ALIVE AGENTS
        
        sarray[0]="Number agents = " + str(pos_count)
        gateway.entry_point.Print_Phrase(sarray[0])

        ## INITIALIZE JVM ARRAYS OF UPDATED AGENT INFORMATION TO PASS TO JAVA
        new_biomass_total = gateway.new_array(gateway.jvm.double,int(pos_count))
        growth_rates_total = gateway.new_array(gateway.jvm.double,int(pos_count))
        index_values_for_java = gateway.new_array(gateway.jvm.int,int(pos_count))
        patch_value_for_java = gateway.new_array(gateway.jvm.double,int(pos_count),num_met_gas+num_met_carbon) # NUMBER OF METABOLITES

        ## FORMAT AGENT INFORMATION FOR MULTIPROCESSING OF AGENT METABOLIC MODELS
        a = np.array((biomass_values,[max_biomass]*pos_count,metabolite_concentrations, index_values, [models]*pos_count, metabolic_states,[rxns]*pos_count,[metabolite_ids]*pos_count,[v_patch]*pos_count,[dt_rxn]*pos_count,[dt_growth]*pos_count),dtype=object).T
        items = list(map(tuple, a))

        ## PERFORM MULTIPROCESSING OPTIMIZATION OF AGENT METABOLIC MODELS
        try:
            results = pool.starmap_async(run_GENRE, items,callback=accumulateResults).get()
        except Exception as e: # RAISE EXCEPTION
            sarray[0]="Multiprocessing FAILED"
            gateway.entry_point.Print_Phrase(sarray[0])
            sarray[0]=str(e)
            gateway.entry_point.Print_Phrase(sarray[0])

        ## FORMAT RESULTS FROM MULTIPROCESSING OPTIMIZATION OF AGENT METABOLIC MODELS
        # CONVERT RESULTS TUPLE TO SEPERATE BIOMASS, METABOLITE CONCENTRATIONS, AND AGENT INDEX ARRAYS
        updated_biomass = [t[0] for t in results]
        updated_biomass = [round(num*1e14,2) for num in updated_biomass]
        updated_metabolite_concentrations = [t[1] for t in results]
        updated_index = [t[2] for t in results]
        growth_rates = [t[3] for t in results]
        growth_rates = [round(num, 10) for num in growth_rates]
        rxn_flux_all = [t[4] for t in results]
        rxn_flux_all = np.array(rxn_flux_all).T
        error = [t[5] for t in results] ## USED TO PRINT ERROR FROM METABOLIC MODEL OPTIMIZATION

        ## FORMAT JVM ARRAYS OF UPDATED AGENT VALUES TO PASS TO JAVA
        for b in range(pos_count):
            new_biomass_total[b] = updated_biomass[b]
            growth_rates_total[b] = growth_rates[b]
            index_values_for_java[b] = int(updated_index[b])
            for m in range((num_met_gas+num_met_carbon)):
                patch_value_for_java[b][m] = updated_metabolite_concentrations[b][m]

        ## SEND AGENT'S UPDATED GROWTH RATE, BIOMASS, AND METABOLITE CONCENTRATIONS TO JAVA
        gateway.entry_point.setGrowthRateFromPython(growth_rates_total, index_values_for_java,int(dead_state),lag_phase,int(t))
        gateway.entry_point.setBiomassFromPython(new_biomass_total, index_values_for_java)
        gateway.entry_point.setPatchFromPython(int(num_met_gas), int(num_met_carbon),patch_value_for_java,index_values_for_java)

        ## RUN METABOLITE DIFFUSION
        gateway.entry_point.Diffuse_Metabolites(int(num_met_gas),initial_gas_concentrations, D_gas,int(num_gas_step),int(num_met_carbon), initial_carbon_concentrations, D_carbon, int(num_carbon_step))
                    
        # SAVE MODEL OUTPUTS: AGENT ATTRIBUTE VALUES AND METABOLITE CONCENTRATION INFORMATION
        if t == 0:
            gateway.Save_cell_info(int(t), int(job_num),output_dir)
        if  t >0:
            gateway.Save_cell_info(int(t), int(job_num),output_dir)
        times_save=[0,25,50,75,100,102,104,106,108,110,112,114,116,118,119] # DEFINE TIMES TO SAVE METABOLITE CONCENTRATIONS
        if  t in times_save:
            gateway.Save_met_info(int(t), int(job_num),output_dir)

        # SAVE MODEL OUTPUTS: AGENT METABOLIC REACTION FLUXES
        df_rxns = pd.DataFrame()
        df_rxns['time'] = np.ones(len(updated_index))*t
        df_rxns['job_num'] = np.ones(len(updated_index))*job_num
        df_rxns['cell index'] = updated_index
        for r,flux in zip(rxns,rxn_flux_all):
            df_rxns[r] = flux

        times_save=[0,25,50,75,100,102,104,106,108,110,112,114,116,118,119] # DEFINE TIMES TO SAVE REACTION FLUXES
        rxns_flux_output_filename = output_dir+'rxns_flux' + str(job_num)+ '.csv'
        if t == 0:
            df_rxns.to_csv(rxns_flux_output_filename)
        if t in times_save:
            df_rxns.to_csv(rxns_flux_output_filename,mode = 'a',header = False)

    ## CLOSE MULTIPROCESSING POOLS
    pool.close();
    pool.join();

## MAIN FUNCTION
if __name__ == '__main__':

    print('MiMICS python file started')

    ################################################ USER INPUTS #######################################################
    ## SET OR IMPORT THE SIMULATION JOB NUMBER
    job_num = 0 # DEFINE JOB NUMBER HERE WHEN RUNNING MIMICS ON PERSONAL COMPUTER
    #job_num=int(os.getenv('NUM_ARRAY')) # COMMENT OUT WHEN RUNNING ON HPC SYSTEM. JOB NUMBER VALUE DEFINED IN MIMICS HPC SLURM JOB FILE.

    ## SET OR IMPORT NUMBER OF MULTIPROCESSING CORES, AND INITIALIZE MULTIPROCESSING WORKERS
    ncpus = 2
    #ncpus=int(os.getenv('NUM_PROCS')) # COMMENT OUT WHEN RUNNING ON HPC SYSTEM. JOB NUMBER VALUE DEFINED IN MIMICS HPC SLURM JOB FILE.

    ## DEFINE FILENAME (XLSX) OF ABM PARAMETER VALUES
    abm_parameters_filename = 'ABM_inputs.xlsx'

    ## DEFINE FILENAME (XLSX) OF NUTRIENT MEDIA
    media_filename = 'Metabolite_inputs.xlsx'

## DEFINE FILENAMES (XML) FOR METABOLIC MODELS
    model0 = "aerobic_state.xml"
    model1 = "denitrification_state.xml"
    model2 = "denitrification_NO_state.xml"
    model3 = "oxidative_stress_state.xml"

    ## DEFINE THE ORDER OF METABOLIC MODELS THAT AGENTS CAN ACCESS VIA THE AGENT'S INTEGER ATTRIBUTE
    metabolic_models_files = [model0,model1,model2,model3]

    ## DEFINE FILENAME (CSV) OF REACTION IDS TO SAVE REACTION FLUX OF EACH AGENT
    rxn_ids_filename = 'MiMICS_short_rxn_list.csv'

    ## DEFINE DIRECTORY FOR SIMULATION OUTPUT FILES
    output_dir = ''
    ############################################# END USER INPUTS ######################################################
    print('User inputs defined')

    ##### IMPORT AND ASSIGN ABM PARAMETERS #####
    abm_parameters = pd.read_excel(abm_parameters_filename, engine='openpyxl')
    abm_parameters = abm_parameters.set_index('Parameter name')['Parameter value'].to_dict() ## SET A DICTIONARY WITH PARAMETER NAMES AND VALUES

    ## SET NUMBER OF SIMULATION TIME STEPS
    total_sim_time = abm_parameters['total_sim_time'] # TOTAL SIMULATION TIME, HOUR
    time_step = abm_parameters['time_step'] # SIMULATION TIME STEP, MINUTES
    num_dt = total_sim_time * 60/time_step; # NUMBER OF SIMULATION TIME STEPS

    ## SET THE X,Y,Z DIMENSIONS OF ABM WORLD AND PDE GRIDS
    patch_length =abm_parameters['patch_scale']
    xdim = abm_parameters['xdim']/patch_length # X-DIMENSION, NUMBER OF PATCHES
    ydim = abm_parameters['ydim']/patch_length  # Y-DIMENSION, NUMBER OF PATCHES
    zdim = abm_parameters['zdim']/patch_length # Z-DIMENSION, NUMBER OF PATCHES

    ## SET ABM AGENT PARAMETERS
    initial_num_agents = abm_parameters['initial_num_agents']
    initial_biomass = abm_parameters['initial_biomass']
    max_biomass = abm_parameters['max_biomass']
    dead_state = abm_parameters['dead_state']
    lag_phase = abm_parameters['lag_phase']
    lag_phase = lag_phase* 60/time_step ## CONVERT LAG PHASE TO SIMULATION TIME STEPS

    ## SET VALUES FOR METABOLIC MODEL OPTIMIZATION
    v_patch = abm_parameters['v_patch'] # VOLUME OF PATCH, L
    dt_rxn = abm_parameters['dt_rxn'] # REACTION TIME SCALE, SECONDS
    dt_growth = abm_parameters['dt_growth'] # BIOMASS GROWTH TIME SCALE, MINUTES

    ## SET METABOLITE REACTION-DIFFUSION TIME SCALES
    dt_diffuse_carbon = abm_parameters['dt_diffuse_carbon']
    dt_diffuse_gas = abm_parameters['dt_diffuse_gas']
    print('ABM parameters imported')

    ##### IMPORT AND DEFINE METABOLITE GRID CONDITIONS #####
    media = pd.read_excel(media_filename, engine='openpyxl')
    metabolite_ids = []
    media_gas = media[media['Gas PDE index'].notna()]
    media_gas = media_gas.sort_values(by='Gas PDE index')
    num_met_gas = len(media_gas)
    metabolite_ids.append(media_gas['Metabolite ID'].values)
    initial_gas_concentrations = media_gas['Metabolite concentration (mM)'].values
    D_gas = media_gas['Diffusion coefficient (cm^2/s)'].values; # GAS DIFFUSION COEFFICIENTS IN BIOFILM

    media_carbon = media[media['Carbon PDE index'].notna()]
    media_carbon = media_carbon.sort_values(by='Carbon PDE index')
    num_met_carbon = len(media_carbon)
    metabolite_ids.append(media_carbon['Metabolite ID'].values)
    initial_carbon_concentrations = media_carbon['Metabolite concentration (mM)'].values
    D_carbon = media_carbon['Diffusion coefficient (cm^2/s)'].values; # CARBON DIFFUSION COEFFICIENTS IN BIOFILM

    metabolite_ids = np.concatenate(metabolite_ids).tolist() # FORMAT METABOLITE IDS

    ## SCALE GAS AND CARBON DIFFUSION COEFFICIENTS
    D_carbon = (D_carbon * dt_diffuse_carbon) / ((patch_length/1e4)**2)
    D_gas = (D_gas* dt_diffuse_gas) / ((patch_length/1e4)**2)

    ## CALCULATE NUMBER OF DIFFUSION TIME STEPS FOR GAS AND CARBON METABOLITES
    num_gas_step = dt_rxn/dt_diffuse_gas
    num_carbon_step = dt_rxn/dt_diffuse_carbon
    print('Metabolite PDE grid conditions imported')

    ##### IMPORT AND DEFINE ORDER OF METABOLIC MODELS #####
    models = []
    for m in metabolic_models_files:
        model = cobra.io.read_sbml_model(m)
        models.append(model)  ## LIST OF METABOLIC MODELS THAT CAN BE USED BY THE AGENTS
    print('Metabolic models imported')

    ##### IMPORT REACTION IDS TO SAVE METABOLIC REACTION FLUXES FOR EACH AGENT #####
    rxns = pd.read_csv(rxn_ids_filename)
    rxns = list(rxns.iloc[:,0])
    print('Reaction IDs imported')

    ############################################### START MiMICS #######################################################
    ## INITIALIZE THE SIMULATION TIME RECORDER
    start_time = time.time()
    
    gateway = JavaGateway(); # INITIALIZE THE PY4J JAVA GATEWAY SERVER

    ## INITIALIZE METABOLITE LIST ARRAYS TO BE PASSED TO JAVA
    initial_gas_concentrations_java = gateway.new_array(gateway.jvm.double,int(num_met_gas)) ## INITIAL CONCENTRATIONS OF GAS METABOLITES
    initial_carbon_concentrations_java = gateway.new_array(gateway.jvm.double,int(num_met_carbon)) ## INITIAL CONCENTRATIONS OF CARBON METABOLITES
    D_gas_java = gateway.new_array(gateway.jvm.double,int(num_met_gas)) ## DIFFUSION COEFFICIENTS OF GAS METABOLITES
    D_carbon_java = gateway.new_array(gateway.jvm.double,int(num_met_carbon)) ## DIFFUSION COEFFICIENTS OF CARBON METABOLITES

    ## SET VALUES IN METABOLITE INITIAL CONCENTRATIONS AND DIFFUSION COEFFICIENTS ARRAYS TO PASS TO JAVA
    for m in range(num_met_gas):
        initial_gas_concentrations_java[m] =initial_gas_concentrations[m]
        D_gas_java[m] =D_gas[m]
    for m in range(num_met_carbon):
        initial_carbon_concentrations_java[m] =initial_carbon_concentrations[m]
        D_carbon_java[m] =D_carbon[m]

    # INITIALIZE ABM AND METABOLITE PDE GRIDS
    gateway.entry_point.run_model0(int(xdim), int(ydim), int(zdim), int(initial_num_agents), initial_biomass, max_biomass,int(num_met_gas), int(num_met_carbon), initial_gas_concentrations_java, initial_carbon_concentrations_java);
    print('Agent-based model and reaction-diffusion model initialized')

    print('Run MiMICS simulation')
    ## RUN MIMICS SIMULATION FOR DESIRED SIMULATION TIME STEPS
    run_MIMICS(ncpus,int(num_dt),models,media,rxns,int(job_num), num_met_gas,initial_gas_concentrations_java, D_gas_java, num_gas_step, num_met_carbon,initial_carbon_concentrations_java, D_carbon_java,num_carbon_step,metabolite_ids,v_patch, dt_rxn, dt_growth,initial_biomass,max_biomass,dead_state,lag_phase,output_dir);
    print('MiMICS simulation finished')

    ## REPORT THE TOTAL SIMULATION TIME
    sarray = gateway.new_array(gateway.jvm.java.lang.String,2) ## STRING JAVA ARRAY TO PRINT STATEMENTS FROM JAVA
    sarray[0]="Total simulation time = " + str(time.time() - start_time) + " seconds"
    gateway.entry_point.Print_Phrase(sarray[0])
    ############################################### END MiMICS #########################################################
