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

## TRACK TIME TO RUN GENRE FOR ALL AGENTS
starttime = time.time();

## FUNCTION TO INITIALIZE MULTIPROCESSING WORKERS
def init():
    try:
        print('Multiprocessing worker initialized', flush=True)
    except Exception: # RAISE EXCEPTION
        sarray = gateway.new_array(gateway.jvm.java.lang.String,2) ## STRING FOR PRINTING FROM JAVA
        sarray[0]="MULTIPROCESSING WORKER FAILED"
        gateway.entry_point.Print_Phrase(sarray[0]) ## PRINT EXCEPTION
        sys.stdout.flush()
        time.sleep(1)

def run_GENRE(old_biomass,oxygen,glucose,index,models, metabolic_state,rxns):

    model = models[metabolic_state] ## ASSIGN METABOLIC MODEL STATE TO AGENT

    ## SET OUTPUT VALUES FOR AGENT (IN CASE EXCEPTION IS RAISED)
    updated_index = index
    new_biomass = old_biomass
    updated_oxygen = oxygen
    updated_glucose = glucose
    growth_rate = 0.0
    rxn_flux_all = np.zeros(len(rxns))

    ## SET BINARY VARIABLES IF MODEL CONTAINS AN EXCHANGE METABOLITE
    model_oxygen_on = 'no'
    model_glucose_on = 'no'
    try:
        vpatch = 2e-16 # PATCH VOLUME, UNITS: L
        dt = 0.1; # SECONDS
        dt = dt/60; # MIN
        dt = dt/60; # HOUR

        if old_biomass > 0:
            conversion_metabolite = vpatch * 1/dt * 1/old_biomass # METABOLITE CONVERSION TERM BETWEEN ABM (mM) > GENRE (mmol/ gdWt hr)

            rxn = "EX_cpd00007_e" ## OXYGEN MODEL SEED ID
            if rxn in model.exchanges:
                model_oxygen_on = 'yes'
                model.reactions.get_by_id(rxn).upper_bound = 1000
                model.reactions.get_by_id(rxn).lower_bound = -oxygen * conversion_metabolite;

            rxn = "EX_cpd00027_e" ## GLUCOSE MODEL SEED ID
            if rxn in model.exchanges:
                model_glucose_on = 'yes'
                model.reactions.get_by_id(rxn).upper_bound = 1000
                model.reactions.get_by_id(rxn).lower_bound = -glucose * conversion_metabolite;

            sol = model.optimize() # OPTIMIZE GENRE

            ## IF SOLUTION IS OPTIMAL
            if  sol.objective_value > 1e-2:
                ## CALCULATE NEW BIOMASS VALUE FOR THIS AGENT
                growth_rate = sol.objective_value # GROWTH RATE (UNITS: 1/hr)
                new_biomass = old_biomass*math.exp(growth_rate*5/60) # UPDATED BIOMASS (UNITS: GRAMS)

                ## CALCULATE NEW METABOLITE CONCENTRATIONS AT THE AGENT'S PATCH
                # CHECK IF ABSOLUTE VALUE OF THE FLUX IS > 1E-12, THEN UPDATE THE METABOLITE CONCENTRATION (UNITS: mM)
                # OXYGEN
                if model_oxygen_on =='yes':
                    ex = 'EX_cpd00007_e'
                    flux = sol.fluxes.loc[ex]
                    if abs(flux) > 1e-12:
                        updated_oxygen = oxygen+flux/conversion_metabolite
                    else:
                        updated_oxygen = oxygen

                    if updated_oxygen < 0:
                        updated_oxygen = 0.0
                # GLUCOSE
                if model_glucose_on =='yes':
                    ex = 'EX_cpd00027_e'
                    flux = sol.fluxes.loc[ex]
                    if abs(flux) > 1e-12:
                        updated_glucose = glucose+flux/conversion_metabolite
                    else:
                        updated_glucose = glucose

                    if updated_glucose < 0:
                        updated_glucose = 0.0
            ## SAVE REACTION FLUX VALUES
            rxn_flux_all = []
            for r in rxns:
                try:
                    rxn_flux = model.reactions.get_by_id(r).flux
                    rxn_flux_all.append(rxn_flux)

                except Exception:
                    rxn_flux = 0.0
                    rxn_flux_all.append(rxn_flux)

    except Exception as e: ## IF GENRE OPTIMIZATION FUNCTION FAILS
        updated_index = index
        new_biomass = old_biomass
        updated_oxygen = oxygen
        updated_glucose = glucose
        growth_rate = 0.0

    return new_biomass, updated_oxygen, updated_glucose, updated_index,growth_rate, tuple(rxn_flux_all)

## FUNCTION TO ACCUMULATE MULTIPROCESSING GENRE RESULTS
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

## FUNCTION TO CALL JAVA FUNCTIONS
def run_MIMICS(num_dt,models,media,rxns,file_name_rxn, job_num,initial_oxygen, D_O2, initial_glucose, D_G):

    sarray = gateway.new_array(gateway.jvm.java.lang.String,2) ## USE TO OUTPUT STATEMENTS FROM PY4J SERVER

    ## IMPORT NUMBER OF MULTIPROCESSING CORES, AND INITIALIZE MULTIPROCESSING WORKERS
    #ncpus=int(os.getenv('NUM_PROCS'))
    ncpus = 2
    pool = Pool(initializer = init, processes=ncpus)

    vpatch = 2e-16 # PATCH VOLUME, UNITS: L
    dt = 0.1; # SECONDS
    dt = dt/60; # MINUTES
    dt = dt/60; # HOURS
    max_biomass = 2e-12;
    conversion_metabolite = vpatch * 1/dt * 1/max_biomass # METABOLITE CONVERSION TERM BETWEEN ABM (mM) > GENRE (mmol/ gdWt hr)

    ## SET LIST OF METABOLIC MODELS THAT CAN BE USED BY THE AGENTS
    models_constrain = [model0, model1]

    # CONSTRAIN METABOLIC MODELS ON NUTRIENT MEDIA (ONLY NEED TO USE FOR FULL iPAU MODEL)
    for model_m in models_constrain:
        medium = {}
        for m,i in zip(media['Metabolite ID'],media.index):
            if 'EX_' + str(m) in model_m.exchanges:
                rxn = 'EX_' + str(m)
                food = media.loc[i]['mM in patch (mmol/L)']
                model_m.reactions.get_by_id(rxn).upper_bound = 1000
                model_m.reactions.get_by_id(rxn).lower_bound = -food*conversion_metabolite
                metrxn = rxn
                medium[metrxn] = food * conversion_metabolite;
        model_m.medium = medium

    ## RUN MULTI-SCALE MODEL FOR DESIRED NUMBER OF SIMULATION TIME STEPS
    for t in range(num_dt):
        sarray = gateway.new_array(gateway.jvm.java.lang.String,2)
        gateway.entry_point.run_model() ## RUN AGENT-BASED MODEL

        sarray = gateway.new_array(gateway.jvm.java.lang.String,2)
        sarray[0]="Time Step = " + str(t)
        gateway.entry_point.Print_Phrase(sarray[0])

        biomass_values = list(gateway.entry_point.getBiomassFromHal()) ## GET AGENT BIOMASS FROM ABM
        biomass_values = [x / 1e14 for x in biomass_values]
        pos_count = len(biomass_values)  ## COUNT NUMBER OF AGENTS
        index_values = list(gateway.entry_point.getIndexFromHal()) ## GET AGENT INDEX FROM ABM
        patch_values_from_java = np.array(gateway.entry_point.getPatchFromHal_All()) ## GET METABOLITE CONCENTRATIONS AT AGENT'S PATCH FROM ABM
        metabolic_states = list(gateway.entry_point.getMetabolicStateFromHal()) ## GET AGENT METABOLIC STATE FROM ABM

        ## FORMAT METABOLITE CONCENTRATIONS
        oxygen = patch_values_from_java[:,0].tolist()
        glucose = patch_values_from_java[:,1].tolist()

        ## INITIALIZE ARRAYS OF UPDATED AGENT INFORMATION TO BE PASSED TO JAVA ABM
        new_biomass_total = gateway.new_array(gateway.jvm.double,int(pos_count))
        growth_rates_total = gateway.new_array(gateway.jvm.double,int(pos_count))
        index_values_for_java = gateway.new_array(gateway.jvm.int,int(pos_count))
        patch_value_for_java = gateway.new_array(gateway.jvm.double,int(pos_count),2) # NUMBER OF METABOLITES

        ## FORMAT AGENT INFORMATION FOR MULTIPROCESSING
        a = np.array((biomass_values,oxygen,glucose, index_values, [models_constrain]*pos_count, metabolic_states,[rxns]*pos_count),dtype=object).T
        items = list(map(tuple, a))
        data_list = list(zip(biomass_values, oxygen, glucose, index_values, [models_constrain]*pos_count,metabolic_states,[rxns]*pos_count))

        ## PERFORM MULTIPROCESSING OPTIMIZATION OF AGENT METABOLIC MODELS
        try:
            results = pool.starmap_async(run_GENRE, data_list,callback=accumulateResults).get()
        except Exception as e:
            sarray[0]="Multiprocessing FAILED"
            gateway.entry_point.Print_Phrase(sarray[0])
            sarray[0]=str(e)
            gateway.entry_point.Print_Phrase(sarray[0])

        ## CONVERT TUPLE TO SEPERATE BIOMASS, METABOLITE CONCENTRATIONS, AND AGENT INDEX ARRAYS
        updated_biomass = [t[0] for t in results]
        updated_biomass = [round(num*1e14,2) for num in updated_biomass]
        updated_oxygen = [t[1] for t in results]
        updated_oxygen = [round(num, 2) for num in updated_oxygen]
        updated_glucose = [t[2] for t in results]
        updated_glucose = [round(num, 2) for num in updated_glucose]
        updated_index = [t[3] for t in results]
        growth_rates = [t[4] for t in results]
        growth_rates = [round(num, 10) for num in growth_rates]
        rxn_flux_all = [t[5] for t in results]
        rxn_flux_all = np.array(rxn_flux_all).T

        ## SAVE AGENT REACTION FLUX VALUES
        df_rxns = pd.DataFrame()
        df_rxns['time'] = np.ones(len(updated_index))*t
        df_rxns['cell index'] = updated_index

        for r,flux in zip(rxns,rxn_flux_all):
            df_rxns[r] = flux

        times_save=[15,30,45,60,75,98,99]
        if t == 0:
            df_rxns.to_csv(file_name_rxn)
        if t in times_save:
            df_rxns.to_csv(file_name_rxn,mode = 'a',header = False)

        ## FORMAT AGENT OUTPUTS FOR JAVA ABM
        for b in range(pos_count):
            new_biomass_total[b] = updated_biomass[b]
            growth_rates_total[b] = growth_rates[b]
            index_values_for_java[b] = int(updated_index[b])
            patch_value_for_java[b][0] = updated_oxygen[b]
            patch_value_for_java[b][1] = updated_glucose[b]

        ## SEND UPDATED GROWTH RATE, BIOMASS, METABOLITE CONCENTRATIONS TO JAVA ABM
        gateway.entry_point.setGrowthRateFromPython(growth_rates_total, index_values_for_java)
        gateway.entry_point.setBiomassFromPython(new_biomass_total, index_values_for_java)
        gateway.entry_point.setPatchFromPython(patch_value_for_java,index_values_for_java)

        # SAVE MODEL OUTPUTS: AGENT AND METABOLITE PATCH INFORMATION
        gateway.entry_point.Diffuse_Metabolites(initial_oxygen, D_O2, initial_glucose, D_G)

        # SAVE MODEL OUTPUTS: AGENT AND METABOLITE PATCH INFORMATION
        if t == 0:
            gateway.Save_cell_info(int(t), int(job_num))
        if  t >0:
            gateway.Save_cell_info(int(t), int(job_num))
        times_save=[0,25,50,91,92,93,94,95,96,97,98,99]
        if  t in times_save:
            gateway.Save_met_info(int(t), int(job_num))

    ## CLOSE MULTIPROCESSING POOLS
    pool.close();
    pool.join();

## MAIN METHOD
if __name__ == '__main__':

    print('MiMICS python file started')

    ## USER INPUT: SET OR IMPORT THE SIMULATION JOB NUMBER
    job_num = 0 # DEFINE HERE WHEN RUNNING MIMICS ON PERSONAL COMPUTER
    #job_num=int(os.getenv('NUM_ARRAY')) # COMMENT OUT WHEN RUNNING ON HPC SYSTEM. JOB NUMBER VALUE DEFINED IN MIMICS HPC SLURM JOB FILE.

    ## USER INPUT: IMPORT METABOLIC MODELS
    model0 = cobra.io.read_sbml_model("insert metabolic model filename")
    model1 = cobra.io.read_sbml_model("insert metabolic model filename")
    models = [model0, model1]  ## SET LIST OF METABOLIC MODELS THAT CAN BE USED BY THE AGENTS

    ## USER INPUT: IMPORT NUTRIENT MEDIA
    media = pd.read_excel('insert nutrient media file', engine='openpyxl')

    print('Metabolic models and nutrient media concentrations imported')

    ## USER INPUT: IMPORT REACTION NAMES TO SAVE REACTION FLUXES FOR EACH AGENT
    rxns = pd.read_csv('insert filename of reaction ids')
    rxns = list(rxns.iloc[:,0])

    ## DEFINE FILE NAME OF REACTION FLUX OUTPUT FILE
    file_name_rxn = 'rxns_flux' + str(job_num)+ '.csv'

    ## USER INPUT: SET NUMBER OF SIMULATION TIME STEPS
    total_sim_time = 0.5 # TOTAL SIMULATION TIME, HOUR
    time_step = 5 # SIMULATION TIME STEP, MINUTES
    num_dt = total_sim_time * 60/time_step; # NUMBER OF SIMULATION TIME STEPS

    ## INITIALIZE THE SIMULATION TIME RECORDER
    start_time = time.time()

    gateway = JavaGateway(); # INITIALIZE THE PY4J JAVA GATEWAY
    sarray = gateway.new_array(gateway.jvm.java.lang.String,2)

    ## USER_INPUT: SET OR IMPORT THE INITIAL NUMBER OF CELLS
    num_cell = 2 # DEFINE HERE WHEN RUNNING MIMICS ON PERSONAL COMPUTER
    #num_cell = int(os.getenv('NUM_CELLS_LIST')) # COMMENT OUT WHEN RUNNING ON HPC SYSTEM. NUM CELLS VALUE DEFINED IN MIMICS HPC SLURM JOB FILE.

    ## USER INPUT: DEFINE INITIAL METABOLITE CONCENTRATIONS
    initial_oxygen = 0.25; # INITIAL OXYGEN CONCENTRATION; UNITS: mM
    initial_glucose = 3.2; # INITIAL GLUCOSE CONCENTRATION; UNITS: mM

    ## USER INPUT: DEFINE METABOLITE DIFFUSION COEFFICIENTS
    D_O2 = 0.5; # SCALED OXYGEN DIFFUSION COEFFICIENT
    D_G = 0.72; # SCALED GLUCOSE DIFFUSION COEFFICIENT

    gateway.entry_point.run_model0(int(num_cell),initial_oxygen, initial_glucose); # INITIALIZE MIMICS MODEL
    print('Agent-based model and reaction-diffusion model initialized')

    print('Run MIMICS simulation')
    run_MIMICS(int(num_dt),models,media,rxns,file_name_rxn,int(job_num), initial_oxygen, D_O2, initial_glucose, D_G);
    print('MIMICS simulation finished')

    ## REPORT THE TOTAL SIMULATION TIME
    sarray[0]="Total simulation time = " + str(time.time() - start_time) + " seconds"
    gateway.entry_point.Print_Phrase(sarray[0])
