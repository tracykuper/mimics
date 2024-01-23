// MIMICS JAVA CLASS: AGENT-BASED MODEL AND REACTION-DIFFUSION MODEL
package MIMICS;

// IMPORT HAL PACKAGES
import HAL.GridsAndAgents.*;
import HAL.Rand;
import HAL.Tools.FileIO;
import HAL.Util;
import static HAL.Util.*;

// IMPORT JAVA PACKAGES
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

// CREATE CELL3D AGENT CLASS
class Cell3D extends SphericalAgent3D<Cell3D, MIMICS>
{
    Random r = new Random();
    Rand rn = new Rand();
    // DEFINE AGENT ATTRIBUTES
    double mass; // BIOMASS
    double forceSum;// AGENT SUM OF FORCES
    double angle; // AGENT ANGLE
    int index; // AGENT INDEX
    double growth_rate; // AGENT GROWTH RATE
    int metabolic_state; // AGENT METABOLIC STATE
    double t_switch_0; // AGENT TRACKING TIME TO SWITCH TO METABOLIC STATE 0
    double t_switch_1; // AGENT TRACKING TIME TO SWITCH TO METABOLIC STATE 1

    // METHODS TO INITIALIZE AN AGENT
    public void Init(double mass,double angle, int state, int index, double t_switch_0, double t_switch_1) {
        // INITIALIZE VALUES OF AGENT ATTRIBUTES
        this.mass = mass;
        this.angle = angle;
        this.radius = G.RADIUS;
        this.index = index;
        this.metabolic_state = state;
        this.t_switch_0 = t_switch_0/ 2;
        this.t_switch_1 = t_switch_1/ 2;
    }

    // METHOD TO CALCULATE SHOVING FORCES BETWEEN AGENTS
    double ForceCalc(double overlap, Cell3D other){
        if(overlap<0) {
            return 0; // IF NO AGENT OVERLAP, RETURN 0
        }
        return G.FORCE_SCALER*overlap; // SCALE OVERLAP (HOOKE'S LAW)
    }

    // METHOD TO PERFORM SHOVING FORCES BETWEEN AGENTS
    public void CalcMove(){
        forceSum=SumForces(G.RADIUS,this::ForceCalc);
    }

    // METHOD TO MOVE AGENTS
    public void MoveDiv(){
        ForceMove();
        ApplyFriction(this.G.FRICTION);
    }

    // METHOD TO SIMULATE AGENT BIOMASS DIVISION
    public int Biomass_Divide(int pop, double initial_biomass, double max_biomass) {
        int birth_counter =0; // VARIABLE TO KEEP TRACK OF NUMBER OF CELL DIVISIONS TO UPDATE AGENT INDEX

        // PERFORM BIOMASS DIVISION FOR AGENTS WITH HIGH BIOMASS
        if (this.mass > max_biomass) {
            Random r = new Random();
            double new_mass = initial_biomass + (max_biomass - initial_biomass)/2 * r.nextDouble(); // CALCULATE DAUGHTER AGENT BIOMASS

            this.mass = this.mass - new_mass; // UPDATE MOTHER AGENT BIOMASS
            double angle2 = this.angle + -10 + (int) (Math.random() * ((10 - (-10)) + 1)); //CALCULATE DAUGHTER AGENT ANGLE

            // CALCULATE DAUGHTER AGENT X,Y,Z COORDINATES
            double dx= 1 + (1 - 0.5) * r.nextDouble();
            double dy= 1 + (1 - 0.5) * r.nextDouble();
            dx= dx * Math.cos(angle);
            dy= dy * Math.sin(angle); ;
            double dz = 0;
            double neg_pos_x = Math.random();
            double neg_pos_y = Math.random();
            if (neg_pos_x > 0.5){
                dx = -dx;
            }
            if (neg_pos_y > 0.5){
                dy = -dy;
            }
            if (Math.random()>0.5 & this.Zpt()<20) {
                Random randomno = new Random();
                dz = randomno.nextGaussian()*0.1+0.6;
                if (dz < 0){
                    dz = 0.4;
                }
            }

            // GENERATE DAUGHTER AGENT
            G.NewAgentPTSafe(this.Xsq()+dx,this.Ysq()+dy,this.Zsq()+dz).Init(new_mass,angle2, this.metabolic_state, pop , this.t_switch_0, this.t_switch_1);
            birth_counter = birth_counter + 1; // UPDATE BIRTH COUNTER
        }
        return birth_counter;
    }

}

// MIMICS MODEL CLASS
public class MIMICS extends AgentGrid3D<Cell3D> {
    public MIMICS(int xdim, int ydim, int zdim) {
        super(xdim, ydim, zdim, Cell3D.class);
        } // CREATE MIMICS CLASS
    Rand rn = new Rand();

    // CREATE ARRAY LISTS FOR AGENT NEIGHBORS
    ArrayList<Cell3D> neighborList = new ArrayList<>();
    ArrayList<Cell3D> neighborList2 = new ArrayList<>();
    ArrayList<double[]> neighborInfo = new ArrayList<>();
    AgentList<Cell3D> cells = new AgentList<>();
    int[] vnHood2d = MooreHood(false); // AGENT 2D NEIGHBORHOOD

    // INITIALIZE MAXIMUM CELL INDEX
    int max_cell_index = 0;

    // INITIALIZE ARRAYS CONTAINING METABOLITE PDE GRIDS
    ArrayList<PDEGrid3D> carbon_metabolites = new ArrayList<PDEGrid3D>(); // CARBON METABOLITE PDE GRIDS
    ArrayList<PDEGrid3D> gas_metabolites = new ArrayList<PDEGrid3D>(); // GAS METABOLITE PDE GRIDS

    // DEFINE AGENT PARAMETERS USED FOR FORCE CALCULATIONS
    double RADIUS; // AGENT RADIUS
    double FORCE_SCALER; // FORCE SCALAR
    double FRICTION; // FRICTION SCALAR

    // METHOD TO SAVE AGENT INFORMATION AT A DESIRED SIMULATION TIME POINT
    public void Save_Cell_Info(int time, int rep,int pop, String file_name) throws IOException {
        double[][] output = new double[pop][11]; // OUTPUT MATRIX: [ROWS: EACH AGENT, COLUMNS: AGENT ATTRIBUTES]
        int c = 0; // INITIALIZE OUTPUT MATRIX ROW INDEX
        for (Cell3D cell : this) { // ITERATE OVER EACH AGENT
            double x = cell.Xpt(); double y = cell.Ypt(); double z = cell.Zpt(); // AGENT X,Y,Z COORDINATES
            // ASSIGN VALUES TO OUTPUT MATRIX
            output[c][0] = time;
            output[c][1] = cell.index;
            output[c][2] = x;
            output[c][3] = y;
            output[c][4] = z;
            output[c][5] = cell.metabolic_state;
            output[c][6] = cell.mass;
            output[c][7] = cell.growth_rate;
            output[c][8] = gas_metabolites.get(0).Get(x, y, z);
            output[c][9] = carbon_metabolites.get(0).Get(x, y, z);
            output[c][10] = rep;
            c = c+1;
        }

        if (time == 0) { // INITIALIZE A NEW CSV OUTPUT FILE IF TIME STEP = 0
        StringBuilder sb = new StringBuilder(); // INITIALIZE STRING BUILDER

        // ASSIGN COLUMN TITLES (correspond to the order of attributes in the output matrix above)
        sb.append("time"); sb.append(",");
        sb.append("cell index"); sb.append(",");
        sb.append("xcor"); sb.append(",");
        sb.append("ycor"); sb.append(",");
        sb.append("zcor"); sb.append(",");
        sb.append("metabolic_state"); sb.append(",");
        sb.append("biomass"); sb.append(",");
        sb.append("growth_rate"); sb.append(",");
        sb.append("oxygen"); sb.append(",");
        sb.append("glucose"); sb.append(",");
        sb.append("job_num"); sb.append('\n');

        // ADD AGENT OUTPUTS TO STRING BUILDER FOR CSV FILE
        for (int j = 0; j < pop; j++) { // ITERATE OVER EACH AGENT
            for (int i = 0; i < 11; i++) { // ITERATE OVER EACH OUTPUT METRIC
                sb.append(output[j][i]); sb.append(","); }
            sb.append('\n'); }

        // WRITE OUTPUT TO NEW CSV
        BufferedWriter br = new BufferedWriter(new FileWriter(file_name));
        br.write(sb.toString());
        br.close();
        }

        if (time != 0) { // APPEND OUTPUT MATRIX TO EXISTING CSV IF TIME IS NOT ZERO
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < pop; j++) { // ITERATE OVER EACH AGENT
                for (int i = 0; i < 11; i++) { // ITERATE OVER EACH OUTPUT METRIC
                    sb.append(output[j][i]); sb.append(","); }
                sb.append('\n');
            }
            // APPEND OUTPUT TO EXISTING CSV
            BufferedWriter br = new BufferedWriter(new FileWriter(file_name,true));
            br.write(sb.toString());
            br.close();
         }
    }

    // METHOD TO SAVE METABOLITE CONCENTRATION INFORMATION
    public void Save_Met_Info(int time, int rep, String file_name) throws IOException {
        double[][] output = new double[xDim * yDim* zDim][7]; // METABOLITE ARRAY: [rows: 3D metabolite grid location, columns: output metrics]
        int c = 0; // INITIALIZE OUTPUT MATRIX ROW INDEX
        // ITERATE OVER EACH PATCH LOCATION IN METABOLITE GRID
        for (int x = 0; x < xDim; x++) {
            for (int z = 0; z < zDim; z++) {
                for (int y = 0; y < yDim; y++) {
                    output[c][0] = x;
                    output[c][1] = y;
                    output[c][2] = z;
                    output[c][3] = gas_metabolites.get(0).Get(x, y, z);
                    output[c][4] = carbon_metabolites.get(0).Get(x, y, z);
                    output[c][5] = time;
                    output[c][6] = rep;
                    c = c + 1;
                }
            }
        }

        if (time == 0) { // INITIALIZE NEW CSV OUTPUT FILE IF TIME STEP = 0
            BufferedWriter br = new BufferedWriter(new FileWriter(file_name));
            StringBuilder sb = new StringBuilder();
            // ASSIGN OUTPUT TITLES (correspond to the order of metrics in output matrix above)
            sb.append("xcor");
            sb.append(",");
            sb.append("ycor");
            sb.append(",");
            sb.append("zcor");
            sb.append(",");
            sb.append("oxygen");
            sb.append(",");
            sb.append("glucose");
            sb.append(",");
            sb.append("time");
            sb.append(",");
            sb.append("job_num"); sb.append('\n');

            // SAVE OUTPUT MATRIX TO STRING BUILDER AND WRITE TO NEW CSV FILE
            for (int j = 0; j < xDim * yDim * zDim; j++) { // ITERATE OVER EACH LOCATION IN 3D METABOLITE GRID
                for (int i = 0; i < 7; i++) { // ITERATE OVER EACH METRIC
                    sb.append(output[j][i]);
                    sb.append(",");
                }
                sb.append('\n');
            }
            br.write(sb.toString());
            br.close();
        }

        if (time != 0) { // IF TIME IS NOT ZERO, APPEND TO EXISTING CSV FILE
            BufferedWriter br = new BufferedWriter(new FileWriter(file_name, true));
            StringBuilder sb = new StringBuilder();

            for (int j = 0; j < xDim * yDim * zDim; j++) { // ITERATE OVER EACH LOCATION IN 3D METABOLITE GRID
                for (int i = 0; i < 7; i++) {
                    sb.append(output[j][i]);
                    sb.append(",");
                }
                sb.append('\n');
            }
            br.write(sb.toString());
            br.close();
        }
    }

    // METHOD TO INITIALIZE METABOLITE CONCENTRATIONS IN ALL LOCATIONS IN METABOLITE PDE GRID
    public void Initialize_Metabolites(int num_met_gas, double[] initial_gas_concentrations, int num_met_carbon,double[] initial_carbon_concentrations) {
        // INITIALIZE GAS METABOLITE CONCENTRATIONS
        for (int m = 0; m < num_met_gas; m++) {
            gas_metabolites.get(m).SetAll(initial_gas_concentrations[m]); // UNITS: mM
            gas_metabolites.get(m).Update();
        }
        // INITIALIZE CARBON METABOLITE CONCENTRATIONS
        for (int m = 0; m < num_met_carbon; m++) {
            carbon_metabolites.get(m).SetAll(initial_carbon_concentrations[m]); // UNITS: mM
            carbon_metabolites.get(m).Update();
        }
    }

    // METHOD TO DIFFUSE GASES
    public void Gas_Diffuse(int num_met_gas,double[] initial_gas_concentration, double[] D_gas, int num_gas_step) {
        for (int i = 0; i < num_gas_step; i++) { // PERFORM GAS DIFFUSION FOR A DISCRETE NUMBER OF TIME STEPS
            // SET AREAS WITHOUT AGENTS TO A CONSTANT CONCENTRATION
            for (int x = 0; x < xDim; x++) {
                for (int z = 0; z < zDim; z++) {
                    for (int y = 0; y < yDim ; y++) {
                        if (GetAgent(x, y, z) == null) {
                            int[] vnHood3d = Util.MooreHood3D(false);
                            int ct = MapEmptyHood(vnHood3d, x, y, z);
                            int ct_max = MapHood(vnHood3d, x, y, z);
                            if (ct == ct_max) {
                                for (int m = 0; m < num_met_gas; m++) {
                                    gas_metabolites.get(m).Set(x, y, z, initial_gas_concentration[m]);
                                }
                            }
                        }
                    }
                }
            }
            for (int m = 0; m < num_met_gas; m++) {
                gas_metabolites.get(m).Update(); // UPDATE METABOLITE CONCENTRATIONS
                gas_metabolites.get(m).DiffusionADI(D_gas[m]); // DIFFUSE METABOLITE CONCENTRATIONS
                gas_metabolites.get(m).Update(); // UPDATE METABOLITE CONCENTRATIONS
            }
        }
    }


    // METHOD TO DIFFUSE CARBON SUBSTRATES
    public void Carbon_Diffuse(int num_met_carbon,double[] initial_carbon_concentration,double[] D_carbon, int num_carbon_step) {
        for (int i = 0; i < num_carbon_step; i++) { // PERFORM CARBON DIFFUSION FOR A DISCRETE NUMBER OF TIME STEPS
            // SET AREAS WITHOUT AGENTS TO A CONSTANT CONCENTRATION
            for (int x = 0; x < xDim; x++) {
                for (int z = 0; z < zDim; z++) {
                    for (int y = 0; y < yDim; y++) {
                        if (GetAgent(x, y, z) == null) {
                            int[] vnHood3d = Util.MooreHood3D(false);
                            int ct = MapEmptyHood(vnHood3d, x, y, z);
                            int ct_max = MapHood(vnHood3d, x, y, z);
                            if (ct == ct_max) {
                                for (int m = 0; m < num_met_carbon; m++) {
                                    carbon_metabolites.get(m).Set(x, y, z, initial_carbon_concentration[m]);
                                }
                            }
                        }
                    }
                }
            }
            for (int m = 0; m < num_met_carbon; m++) {
                carbon_metabolites.get(m).Update(); // UPDATE METABOLITE CONCENTRATIONS
                carbon_metabolites.get(m).DiffusionADI(D_carbon[m]); // DIFFUSE METABOLITE CONCENTRATIONS
                carbon_metabolites.get(m).Update(); // UPDATE METABOLITE CONCENTRATIONS
            }
        }
    }

    // METHOD TO PERFORM AGENT METHODS
    public void StepCells(double initial_biomass,double max_biomass,int dead_state) {
        CleanAgents(); // CLEAN REMOVED AGENTS IF NECESSARY
        // ITERATE OVER EACH AGENT
        for (Cell3D cell : this) {
            if (cell.metabolic_state != dead_state) { // IF THE CELL IS NOT DEAD

                // PERFORM AGENT BIOMASS DIVISION
                int birth_counter = cell.Biomass_Divide(max_cell_index,initial_biomass,max_biomass);

                //  RESET MAXIMUM CELL INDEX (USED TO ASSIGN NEW CELL INDEX TO DAUGHTER CELLS) IF DIVISION EVENT OCCURED
                if (birth_counter > 0) {
                    max_cell_index = max_cell_index + 1;
                }

                // PERFORM CELL MECHANICS
                for (int i = 0; i < 5; i++) {
                    cell.CalcMove();
                    cell.MoveDiv();
                }
            }
        }
            CleanAgents(); // CLEAN REMOVED AGENTS IF NECESSARY
    }

    // METHOD TO RANDOMLY INITIALIZE AGENTS
    public void Initialize_Random(int b0, double initial_biomass, double max_biomass) {
        max_cell_index = b0; // SET MAXIMUM CELL INDEX
        for (int j = 0; j <= max_cell_index-1; j++) { // j DEFINES THE AGENT INDEX
            double x = rn.Double(xDim); // X-COORDINATE
            double y = rn.Double(yDim); // Y-COORDINATE
            Random r = new Random();
            double old_mass =  initial_biomass + (max_biomass - initial_biomass) * r.nextDouble(); // SET AGENT BIOMASS
            int angle = rn.Int(360); // INITIALIZE CELL ANGLE
            int metabolic_state = 0; // INITIALIZE METABOLIC MODEL STATE
            NewAgentPT(x, y, 0).Init(old_mass,angle, metabolic_state, j, 0, 0);
        }
    }

    // METHOD TO CALL MAIN METHOD
    static int count = 0;
    public static void mainCaller() throws IOException, NoSuchFieldException, IllegalAccessException {
        count++;
        while (count == 1) {
            main(null);
        }
    }

    public static void main(String[] args) throws IOException, NoSuchFieldException, IllegalAccessException {
        mainCaller(); // MAIN METHOD USED EXTERNALLY IN MIMICS GATEWAY CLASS
    }
}