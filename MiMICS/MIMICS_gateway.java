// MIMICS PY4J GATEWAY FILE
package MIMICS;

//IMPORT JAVA PACKAGES
import HAL.GridsAndAgents.PDEGrid3D;
import py4j.GatewayServer; //PY4J PACKAGE
import java.io.*;

// CREATE MIMICS CLASS (NAMED BIS)
public class MIMICS_gateway {
    public static MIMICS bis;
    public MIMICS_gateway() {
    }

    // FUNCTION TO TO RETURN MIMICS CLASS (EXECUTED FROM PYTHON)
    public MIMICS getStack() {
        return bis;
    }

    // FUNCTION TO PRINT PHRASE (EXECUTED FROM PYTHON)
    public static void Print_Phrase(String phrase) {
        System.out.println(phrase);
    }

    // FUNCTION TO INITIALIZE MIMICS (EXECUTED FROM PYTHON)
    public static void run_model0(int xdim, int ydim, int zdim, int b0, double initial_biomass,double max_biomass, int num_met_gas, int num_met_carbon, double[] initial_gas_concentrations, double[] initial_carbon_concentrations) {
        bis = new MIMICS(xdim, ydim, zdim); //CREATE MIMICS OBJECT AND SET X,Y,Z DIMENSIONS
        bis.wrapX = true; // SET PERIODIC BOUNDARY CONDITIONS
        bis.wrapY = true; // SET PERIODIC BOUNDARY CONDITIONS

        // INITIALIZE RANDOM DISTRIBUTION OF AGENTS IN ABM
        bis.Initialize_Random(b0, initial_biomass,max_biomass);

        // INITIALIZE METABOLITE PDE GRIDS
        for (int m = 0; m< num_met_gas; m++) { //GASEOUS METABOLITES
            bis.gas_metabolites.add(new PDEGrid3D(xdim, ydim, zdim));
        }
        for (int m = 0; m< num_met_carbon; m++) { //CARBON METABOLITES
            bis.carbon_metabolites.add(new PDEGrid3D(xdim, ydim, zdim));
        }

        // INITIALIZE METABOLITE CONCENTRATIONS
        bis.Initialize_Metabolites(num_met_gas,initial_gas_concentrations,num_met_carbon,initial_carbon_concentrations);
    }

    // DIFFUSE METABOLITES
    public static void Diffuse_Metabolites(int num_met_gas,double[] initial_gas_concentrations, double[] D_gas, int num_gas_step,int num_met_carbon,double[] initial_carbon_concentrations, double[] D_carbon, int num_carbon_step) {
        bis.Gas_Diffuse(num_met_gas,initial_gas_concentrations, D_gas, num_gas_step); // GAS METABOLITE DIFFUSION
        bis.Carbon_Diffuse(num_met_carbon, initial_carbon_concentrations,D_carbon, num_carbon_step); // CARBON METABOLITE DIFFUSION
    }

    // FUNCTION TO RUN AGENT-BASED MODEL (EXECUTED FROM PYTHON)
    public static void run_ABM(double initial_biomass, double max_biomass,int dead_state) {
        bis.StepCells(initial_biomass,max_biomass,dead_state);
    }

    // FUNCTION TO SAVE METABOLITE OUTPUTS (EXECUTED FROM PYTHON)
    public static void Save_met_info(int time, int rep,String dir) throws IOException {
        bis.Save_Met_Info(time,rep,dir + "met_grid" + Integer.toString(rep) + ".csv"); // SAVE METABOLITE PATCH INFORMATION
    }

    // FUNCTION TO SAVE AGENT OUTPUTS (EXECUTED FROM PYTHON)
    public static void Save_cell_info(int time, int rep,String dir) throws IOException {
        bis.Save_Cell_Info(time,rep,bis.Pop(),dir + "agent_properties" + Integer.toString(rep) + ".csv") ; // SAVE AGENT INFORMATION
    }

    // FUNCTION TO SEND AGENT BIOMASS FROM ABM TO PYTHON (EXECUTED FROM PYTHON)
    public double[] getBiomassFromHal() {
        double[] biomassFromHal = new double[bis.max_cell_index]; // AGENT BIOMASS ARRAY
        // ITERATE OVER ALL AGENTS
        for (Cell3D cell : bis.AllAgents()) {
            biomassFromHal[cell.index] = Math.round(cell.mass * 1e14 *1000)/1000; // BIOMASS ASSIGNED TO RESPECTIVE AGENT INDEX IN ARRAY
        }
        return biomassFromHal;
    }

    // FUNCTION TO SEND AGENT INDEX FROM ABM TO PYTHON (EXECUTED FROM PYTHON)
    public int[] getIndexFromHal() {
        int[] indexFromHal = new int[bis.max_cell_index]; // AGENT INDEX ARRAY
        // ITERATE OVER ALL AGENTS
        for (Cell3D cell : bis.AllAgents()) {
            indexFromHal[cell.index] = cell.index; // INDEX ASSIGNED TO RESPECTIVE AGENT INDEX IN ARRAY
        }
        return indexFromHal;
    }

    // FUNCTION TO SEND AGENT METABOLIC STATE ATTRIBUTE FROM ABM TO PYTHON (EXECUTED FROM PYTHON)
    public int[] getMetabolicStateFromHal() {
        int[] metabolicstatesFromHal = new int[bis.max_cell_index]; // AGENT METABOLIC STATE ARRAY
        // ITERATE OVER ALL AGENTS
        for (Cell3D cell : bis.AllAgents()) {
            metabolicstatesFromHal[cell.index]=cell.metabolic_state; // METABOLIC STATE ASSIGNED TO RESPECTIVE AGENT INDEX IN ARRAY
        }
        return metabolicstatesFromHal;
    }

    // FUNCTION TO SEND METABOLITE CONCENTRATIONS AT THE AGENT'S LOCATION FROM ABM TO PYTHON (EXECUTED FROM PYTHON)
    public double[][] getPatchFromHal_All(int num_met_gas, int num_met_carbon) {
        double[][] patchFromHal = new double[bis.max_cell_index][num_met_gas + num_met_carbon]; // AGENT METABOLITE CONCENTRATION ARRAY, EACH ROW CORRESPONDS TO AN AGENT INDEX, EACH COLUMN CORRESPONDS TO A METABOLITE
        // ITERATE OVER ALL AGENTS AND GET METABOLITE CONCENTRATIONS OF THE AGENT'S PATCH
        // ROUND METABOLITES TO CONSERVE DATA PASSED BETWEEN JAVA AND PYTHON
        for (Cell3D cell : bis.AllAgents()) {
            //GAS METABOLITES
            for (int m = 0; m < num_met_gas; m++) {
                double gas_here = (double) Math.round(bis.gas_metabolites.get(m).Get(cell.Xsq(), cell.Ysq(), cell.Zsq()) * 1000) / 1000;
                patchFromHal[cell.index][m] = gas_here;
            }
            //CARBON METABOLITES
            for (int m = 0; m < num_met_carbon; m++) {
                double carbon_here = (double) Math.round(bis.carbon_metabolites.get(m).Get(cell.Xsq(), cell.Ysq(), cell.Zsq()) * 1000) / 1000;
                patchFromHal[cell.index][m+num_met_gas] = carbon_here;
            }
        }
        return patchFromHal;
    }

    // FUNCTION TO UPDATE METABOLITE CONCENTRATIONS AT AGENT'S LOCATIONS; PASSED FROM PYTHON TO ABM (EXECUTED FROM PYTHON)
    public void setPatchFromPython(int num_met_gas, int num_met_carbon, double[][] met, int[] index) {
        // ITERATE OVER AGENT INDICES PASSED FROM PYTHON
        for (int g = 0; g <index.length; g++) {
            // ITERATE OVER AGENTS IN ABM
            for (Cell3D cell : bis.AllAgents()) {
                 if (cell.index ==index[g]) { // CHECK IF AGENT INDEX IS EQUAL TO INDEX FROM ARRAY IN PYTHON
                     // SET GAS METABOLITE CONCENTRATIONS AT AGENT'S LOCATION
                     for (int m = 0; m < num_met_gas; m++) {
                         cell.G.gas_metabolites.get(m).Set(cell.Xsq(), cell.Ysq(), cell.Zsq(), met[g][m]);
                     }
                     // SET CARBON METABOLITE CONCENTRATIONS AT AGENT'S LOCATION
                     for (int m = 0; m < num_met_carbon; m++) {
                         cell.G.carbon_metabolites.get(m).Set(cell.Xsq(), cell.Ysq(), cell.Zsq(), met[g][m+num_met_gas]);
                     }
                 }
            }
        }
        // UPDATE CONCENTRATIONS IN GAS AND CARBON METABOLITE GRIDS
        for (int m = 0; m < num_met_carbon; m++) {
            bis.carbon_metabolites.get(m).Update();
        }
        for (int m = 0; m < num_met_gas; m++) {
            bis.gas_metabolites.get(m).Update();
        }
    }

    // FUNCTION TO UPDATE AGENT'S GROWTH RATE; PASSED FROM PYTHON TO ABM (EXECUTED FROM PYTHON)
    public void setGrowthRateFromPython(double[] growth_rate, int[] index, int dead_state) {
        // ITERATE OVER AGENT INDICES
        for (int g = 0; g < index.length; g++) {
            // ITERATE OVER AGENTS IN ABM
            for (Cell3D cell : bis.AllAgents()) {
                if (cell.index == index[g]) { // CHECK IF AGENT INDEX IS EQUAL TO INDEX FROM PYTHON
                    cell.growth_rate = growth_rate[g]; // UPDATE AGENT BIOMASS

                    // ASSIGN DEAD CELL STATE
                    if (cell.growth_rate ==0){
                        cell.metabolic_state = dead_state;
                    }
                }
            }
        }
    }

    // FUNCTION TO UPDATE AGENT'S GROWTH RATE FROM PYTHON TO ABM (EXECUTED FROM PYTHON)
    public void setBiomassFromPython(double[] mass, int[] index) {
        // ITERATE OVER AGENT INDICES
        for (int g = 0; g < index.length; g++) {
            // ITERATE OVER AGENTS IN ABM
            for (Cell3D cell : bis.AllAgents()) { // CHECK IF AGENT INDEX IS EQUAL TO INDEX FROM PYTHON
                if (cell.index == index[g]) {
                    cell.mass = mass[g] / 1e14; // DIVIDE BY 1e14 (used to reduce bytes)
                }
            }
        }
    }

    // MIMICS GATEWAY MAIN METHOD
    public static void main(String[] args) throws IOException, InterruptedException {
        GatewayServer gatewayServer = new GatewayServer(new MIMICS_gateway()); // INITIALIZE JAVA PY4J GATEWAY SERVER
        gatewayServer.start(); // START JAVA GATEWAY SERVER FOR PYTHON ACCESS
        System.out.println("PY4J Gateway Server Started");
//
//        // USE  CODE BELOW FOR MIMICS GATEWAY SERVER TO AUTOMATICALLY RUN MIMICS PYTHON FILE (agents_FBA.py)
//
//        String command = "python agents_FBA.py"; // PYTHON FILE TO OPTIMIZE EACH AGENT'S GENRE
//        Process p = Runtime.getRuntime().exec(command); // CALL PYTHON FILE
//        p.waitFor();
//
//        // PRINT STATEMENTS CALLED FROM THE PYTHON FILE
//        try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
//            String line;
//            while ((line = br.readLine()) != null) {
//                System.out.println(line);
//            }
//        }
//        gatewayServer.shutdown(); //SHUT DOWN JAVA GATEWAY SERVER WHEN MIMICS IS FINISHED
    }
}

