// JAVA GATEWAY FILE

package MIMICS;

//IMPORT JAVA PACKAGES
import py4j.GatewayServer; //PY4J PACKAGE
import java.io.*;

// CREATE AGENT-BASED MODEL CLASS (BIS)
public class MIMICS_gateway {
    public static MIMICS bis;
    public MIMICS_gateway(int xdim, int ydim, int zdim) {
        bis = new MIMICS(xdim, ydim, zdim); //CREATE ABM OBJECT AND SET DIMENSIONS
        bis.wrapX = true; // SET PERIODIC BOUNDARY CONDITIONS
        bis.wrapY = true; // SET PERIODIC BOUNDARY CONDITIONS
    }

    // FUNCTION TO TO RETURN AGENT-BASED MODEL (EXECUTED FROM PYTHON)
    public MIMICS getStack() {
        return bis;
    }

    // FUNCTION TO PRINT PHRASE (EXECUTED FROM PYTHON)
    public static void Print_Phrase(String phrase) {
        System.out.println(phrase);
    }

    // FUNCTION TO INITIALIZE MODEL (EXECUTED FROM PYTHON)
    // b0: user-input of initial number of agents, or the microscopy image replicate
    public static void run_model0(int b0, double initial_oxygen, double initial_glucose) {
        bis.ResetHard(); //RESET MODEL

        // RANDOM DISTRIBUTION OF INITIAL AGENTS
        bis.Initialize_Random(b0);

        // INITIALIZE METABOLITE CONCENTRATIONS IN ABM PATCHES
        bis.Initialize_Metabolites(initial_oxygen,initial_glucose);
    }

    public static void Diffuse_Metabolites(double initial_oxygen, double D_O2, double initial_glucose, double D_G) {
        // DIFFUSE METABOLITES
        bis.Gas_Diffuse(initial_oxygen, D_O2);
        bis.Carbon_Diffuse(initial_glucose,D_G);
    }

    // FUNCTION TO RUN AGENT-BASED MODEL (EXECUTED FROM PYTHON)
    public static void run_model() {
        bis.StepCells();
    }

    // FUNCTION TO SAVE MODEL OUTPUTS (EXECUTED FROM PYTHON)
    public static void Save_met_info(int time, int rep) throws IOException {
        String file_name = "met_grid"+ Integer.toString((rep))+".csv";
        bis.Save_Met_Info(time,rep,file_name); // SAVE METABOLITE PATCH INFORMATION
    }

    public static void Save_cell_info(int time, int rep) throws IOException {
        String file_name ="agent_properties"+ Integer.toString((rep))+".csv";
        bis.Save_Cell_Info(time,rep,bis.Pop(),file_name); // SAVE AGENT INFORMATION
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

    // FUNCTION TO SEND AGENT PHENOTYPE FROM ABM TO PYTHON (EXECUTED FROM PYTHON)
    public int[] getMetabolicStateFromHal() {
        int[] metabolicstatesFromHal = new int[bis.max_cell_index]; // AGENT PHENOTYPE ARRAY
        // ITERATE OVER ALL AGENTS
        for (Cell3D cell : bis.AllAgents()) {
            metabolicstatesFromHal[cell.index]=cell.metabolic_state; // PHENOTYPE ASSIGNED TO RESPECTIVE AGENT INDEX IN ARRAY
        }
        return metabolicstatesFromHal;
    }


    // FUNCTION TO SEND METABOLITE CONCENTRATIONS OF THE AGENT'S PATCH FROM ABM TO PYTHON (EXECUTED FROM PYTHON)
    public double[][] getPatchFromHal_All() {
        double[][] patchFromHal = new double[bis.max_cell_index][3]; // AGENT METABOLITE CONCENTRATION ARRAY, EACH ROW CORRESPONDS TO AN AGENT, EACH COLUMN CORRESPONDS TO A METABOLITE
        // ITERATE OVER ALL AGENTS AND GET METABOLITE CONCENTRATIONS OF THE AGENT'S PATCH
        // ROUND METABOLITES TO CONSERVE DATA PASSED BETWEEN JAVA AND PYTHON
        for (Cell3D cell : bis.AllAgents()) {
            double o2_here = (double) Math.round( bis.oxygen.Get(cell.Xsq(), cell.Ysq(), cell.Zsq()) *1000)/1000;
            patchFromHal[cell.index][0] = o2_here;

            double glucose_here = (double) Math.round( bis.glucose.Get(cell.Xsq(), cell.Ysq(), cell.Zsq()) *1000)/1000;
            patchFromHal[cell.index][1] = glucose_here;

        }
        return patchFromHal;
    }



    // FUNCTION TO UPDATE METABOLITE CONCENTRATIONS OF AGENT'S PATCH FROM PYTHON TO ABM (EXECUTED FROM PYTHON)
    public void setPatchFromPython(double[][] met, int[] index) {
        // ITERATE OVER AGENT INDICES
        for (int g = 0; g <index.length; g++) {
            // ITERATE OVER AGENTS IN ABM
            for (Cell3D cell : bis.AllAgents()) {
                 if (cell.index ==index[g]) { // CHECK IF AGENT INDEX IS EQUAL TO INDEX FROM PYTHON
                     // SET METABOLITE CONCENTRATIONS IN AGENT PATCH
                    cell.G.oxygen.Set(cell.Xsq(), cell.Ysq(), cell.Zsq(), met[g][0]);
                    cell.G.glucose.Set(cell.Xsq(), cell.Ysq(), cell.Zsq(), met[g][1]);
                 }
            }
        }
        // UPDATE METABOLITE GRIDS
        bis.oxygen.Update();
        bis.glucose.Update();
    }


    // FUNCTION TO UPDATE AGENT'S GROWTH RATE FROM PYTHON TO ABM (EXECUTED FROM PYTHON)
    public void setGrowthRateFromPython(double[] growth_rate, int[] index) {
        // ITERATE OVER AGENT INDICES
        for (int g = 0; g < index.length; g++) {
            // ITERATE OVER AGENTS IN ABM
            for (Cell3D cell : bis.AllAgents()) {
                if (cell.index == index[g]) { // CHECK IF AGENT INDEX IS EQUAL TO INDEX FROM PYTHON
                    cell.growth_rate = growth_rate[g]; // UPDATE AGENT BIOMASS
                    if (cell.growth_rate ==0){
                        cell.metabolic_state = 10;
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
                    cell.mass = mass[g] / 1e14; // DIVIDE BY 1e13 (used to reduce bytes)
                }
            }
        }
    }

    // MAIN METHOD
    public static void main(String[] args) throws IOException, InterruptedException {
        // USER INPUT: DEFINE 3D WORLD DIMENSIONS
        int xdim = (int) 50; // X-DIMENSION, UNITS: MICROMETERS
        int ydim = (int) 50; // Y-DIMENSION, UNITS: MICROMETERS
        int zdim = (int) 10; // Z-DIMENSION, UNITS: MICROMETERS

        GatewayServer gatewayServer = new GatewayServer(new MIMICS_gateway(xdim,ydim,zdim)); // INITIALIZE JAVA GATEWAY SERVER, CALLING ABM JAVA CLASS
        gatewayServer.start(); // START JAVA GATEWAY SERVER FOR PYTHON ACCESS
        System.out.println("PY4J Gateway Server Started");

        // USE BELOW CODE TO AUTOMICALLY RUN PYTHON FILE WHEN MIMICS GATWAY IS STARTED
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
//        gatewayServer.shutdown(); //SHUT DOWN JAVA GATEWAY SERVER WHEN MULTI-SCALE MODEL IS FINISHED
    }
}

