// MiMICS: AGENT-BASED MODEL AND REACTION-DIFFUSION MODEL
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

// CREATE AGENT CLASS
class Cell3D extends SphericalAgent3D<Cell3D, MIMICS>
{
    Random r = new Random();
    Rand rn = new Rand();
    // DEFINE AGENT VARIABLES
    double mass; // BIOMASS
    double forceSum;// AGENT SUM OF FORCES
    double angle; // AGENT ANGLE
    int index; // AGENT INDEX
    double growth_rate; // AGENT GROWTH RATE
    int metabolic_state; // AGENT METABOLIC STATE
    double t_switch_0; // AGENT TRACKING TIME TO SWITCH TO METABOLIC STATE 0
    double t_switch_1; // AGENT TRACKING TIME TO SWITCH TO METABOLIC STATE 1
    int pili_label;

    // FUNCTION TO INITIALIZE AN AGENT
    public void Init(double mass, double angle, int state, int index, double t_switch_0, double t_switch_1) {
        // INITIALIZE AGENT PROPERTIES
        this.mass = mass;
        this.angle = angle;
        this.radius = G.RADIUS;
        this.index = index;
        this.metabolic_state = state;
        this.t_switch_0 = t_switch_0/ 2;
        this.t_switch_1 = t_switch_1/ 2;
    }

    // FUNCTION TO CALCULATE SHOVING FORCES BETWEEN AGENTS
    double ForceCalc(double overlap, Cell3D other){
        if(overlap<0) {
            return 0; // IF NO AGENT OVERLAP, RETURN 0
        }
        return G.FORCE_SCALER*overlap; // SCALE OVERLAP (HOOKE'S LAW)
    }

    // FUNCTION TO PERFORM SHOVING FORCES BETWEEN AGENTS
    public void CalcMove(){
        forceSum=SumForces(G.RADIUS,this::ForceCalc);
    }

    // FUNCTION TO MOVE AGENTS
    public void MoveDiv(){
        ForceMove();
        ApplyFriction(this.G.FRICTION);
    }

    // FUNCTION TO SIMULATE AGENT BIOMASS GROWTH AND DIVISION
    public int Biomass_Divide(int pop) {
        int birth_counter =0; // KEEP TRACK OF NUMBER OF CELL DIVISIONS TO UPDATE AGENT INDEX

        // PERFORM AGENT DIVISION FOR AGENTS WITH HIGH BIOMASS
        double max_biomass =2 * 1e-12;
        if (this.mass > max_biomass) {
            Random r = new Random();
            double new_mass = 1e-12 + (0.5e-12) * r.nextDouble();

            this.mass = this.mass - new_mass; //UPDATE MOTHER AGENT BIOMASS
            double angle2 = this.angle + -10 + (int) (Math.random() * ((10 - (-10)) + 1)); //CALCULATE DAUGHTER AGENT ANGLE

            //CALCULATE DAUGHTER AGENT X,Y,Z COORDINATES
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

            //GENERATE DAUGHTER AGENT
            G.NewAgentPTSafe(this.Xsq()+dx,this.Ysq()+dy,this.Zsq()+dz).Init(new_mass, angle2, this.metabolic_state, pop , this.t_switch_0, this.t_switch_1);
            birth_counter = birth_counter + 1; //UPDATE BIRTH COUNTER

        }
        return birth_counter;
    }


    // FUNCTION TO CHANGE AGENT METABOLIC STATE BASED ON MECHANISTIC RULES
    public void Change_metabolic_state() {
            }
}

// HAL MODEL CLASS
public class MIMICS extends AgentGrid3D<Cell3D> {
    public MIMICS(int xdim, int ydim, int zdim) {
        super(xdim, ydim, zdim, Cell3D.class); } // CREATE HAL MODEL CLASS
    Rand rn = new Rand();

    // CREATE ARRAY LISTS FOR AGENT NEIGHBORS
    ArrayList<Cell3D> neighborList = new ArrayList<>();
    ArrayList<Cell3D> neighborList2 = new ArrayList<>();
    ArrayList<double[]> neighborInfo = new ArrayList<>();
    AgentList<Cell3D> cells = new AgentList<>();
    int[] vnHood2d = MooreHood(false); //AGENT 2D NEIGHBORHOOD

    // INITIALIZE MAXIMUM CELL INDEX
    int max_cell_index =0;

    // USER INPUT: INITIALIZE PDE GRIDS FOR METABOLITES IN BIOFILM
    PDEGrid3D oxygen = new PDEGrid3D(xDim, yDim,zDim);
    PDEGrid3D glucose = new PDEGrid3D(xDim, yDim,zDim);

    // DEFINE AGENT PARAMETERS USED FOR FORCE CALCULATIONS
    double RADIUS = 0.45; // AGENT RADIUS
    double FORCE_SCALER = 1; // FORCE SCALAR
    double FRICTION = 0.5; // FRICTION SCALAR, 1 = no friction, 0 = maximum friction and x,y velocity will be 0

    // FUNCTION TO SAVE AGENT INFORMATION AT A DESIRED TIME POINT
    public void Save_Cell_Info(int time, int rep,int pop, String file_name) throws IOException {
        double[][] output = new double[pop][13]; // OUTPUT MATRIX: [ROWS: EACH AGENT, COLUMNS: AGENT PROPERTIES]
        int c = 0; // INITIALIZE OUTPUT MATRIX INDEX
        for (Cell3D cell : this) { // ITERATE OVER EACH AGENT
            double x = cell.Xpt(); double y = cell.Ypt(); double z = cell.Zpt(); // AGENT COORDINATES
            output[c][0] = time;
            output[c][1] = cell.index;
            output[c][2] = x;
            output[c][3] = y;
            output[c][4] = z;
            output[c][5] = cell.metabolic_state;
            output[c][6] = cell.mass;
            output[c][7] = cell.growth_rate;
            output[c][8] = cell.pili_label;
            output[c][9] = oxygen.Get(x, y, z);
            output[c][10] = glucose.Get(x, y, z);
            output[c][12] = rep;
            c = c+1;
        }

        if (time == 0) { // INITIALIZE NEW CSV OUTPUT FILE IF TIME STEP = 0
        StringBuilder sb = new StringBuilder(); // INITIALIZE STRING BUILDER
        // ASSIGN OUTPUT TITLES (correspond to order of metrics in output)
        sb.append("time"); sb.append(",");
        sb.append("cell index"); sb.append(",");
        sb.append("xcor"); sb.append(",");
        sb.append("ycor"); sb.append(",");
        sb.append("zcor"); sb.append(",");
        sb.append("metabolic_state"); sb.append(",");
        sb.append("biomass"); sb.append(",");
        sb.append("growth_rate"); sb.append(",");
        sb.append("pili_label"); sb.append(",");
        sb.append("oxygen"); sb.append(",");
        sb.append("glucose"); sb.append(",");
        sb.append("FOV"); sb.append('\n');

        // ADD AGENT OUTPUTS TO STRING BUILDER FOR CSV FILE
        for (int j = 0; j < pop; j++) { // ITERATE OVER EACH AGENT
            for (int i = 0; i < 13; i++) { // ITERATE OVER EACH OUTPUT METRIC
                sb.append(output[j][i]); sb.append(","); }
            sb.append('\n'); }

        // WRITE OUTPUT TO NEW CSV
        BufferedWriter br = new BufferedWriter(new FileWriter(file_name));
        br.write(sb.toString());
        br.close();
        }

        if (time != 0) { // APPEND TO EXISTING CSV IF TIME IS NOT ZERO
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < pop; j++) { // ITERATE OVER EACH AGENT
                for (int i = 0; i < 13; i++) { // ITERATE OVER EACH OUTPUT METRIC
                    sb.append(output[j][i]); sb.append(","); }
                sb.append('\n');
            }
            // APPEND OUTPUT TO EXISTING CSV
            BufferedWriter br = new BufferedWriter(new FileWriter(file_name,true));
            br.write(sb.toString());
            br.close();
         }
    }

    // FUNCTION TO SAVE METABOLITE PATCH INFORMATION
    public void Save_Met_Info(int time, int rep, String file_name) throws IOException {
        double[][] output = new double[xDim * yDim* zDim][6]; // METABOLITE ARRAY: [rows: 3D metabolite grid location, columns: output metrics]
        int c = 0; // OUTPUT ROW INDEX
        // ITERATE OVER EACH PATCH LOCATION IN METABOLITE GRID
        for (int x = 0; x < xDim; x++) {
            for (int z = 0; z < zDim; z++) {
                for (int y = 0; y < yDim; y++) {
                    output[c][0] = x;
                    output[c][1] = y;
                    output[c][2] = z;
                    output[c][3] = oxygen.Get(x, y, z);
                    output[c][4] = glucose.Get(x, y, z);
                    output[c][5] = time;
                    c = c + 1;
                }
            }
        }

        if (time == 0) { // INITIALIZE NEW CSV OUTPUT FILE IF TIME STEP = 0
            BufferedWriter br = new BufferedWriter(new FileWriter(file_name));
            StringBuilder sb = new StringBuilder();
            // ASSIGN OUTPUT TITLES (correspond to order of metrics in output matrix)
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
            sb.append('\n');

            // SAVE OUTPUT MATRIX TO STRING BUILDER AND WRITE TO NEW CSV FILE
            for (int j = 0; j < xDim * yDim * zDim; j++) { // ITERATE OVER EACH PATCH LOCATION IN 3D METABOLITE GRID
                for (int i = 0; i < 6; i++) { // ITERATE OVER EACH METRIC
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

            for (int j = 0; j < xDim * yDim * zDim; j++) { // ITERATE OVER EACH PATCH LOCATION IN 3D METABOLITE GRID
                for (int i = 0; i < 6; i++) {
                    sb.append(output[j][i]);
                    sb.append(",");
                }
                sb.append('\n');
            }
            br.write(sb.toString());
            br.close();
        }
    }

    // FUNCTION TO INITIALIZE METABOLITE CONCENTRATIONS IN ALL PATCHES
    public void Initialize_Metabolites(double initial_oxygen, double initial_glucose) {
        // USER INPUT: DEFINE INITIAL METABOLITE CONCENTRATIONS
        oxygen.SetAll(initial_oxygen); // UNITS: mM
        oxygen.Update();

        glucose.SetAll(initial_glucose); // UNITS: mM
        glucose.Update();
    }

    // FUNCTION TO DIFFUSE GASES
    public void Gas_Diffuse(double initial_oxygen, double D_O2) {
        for (int i = 0; i < 20; i++) {
            for (int x = 0; x < xDim; x++) {
                for (int z = 0; z < zDim; z++) {
                    for (int y = 0; y < yDim ; y++) {
                        if (GetAgent(x, y, z) == null) {
                            int[] vnHood3d = Util.MooreHood3D(false);
                            int ct = MapEmptyHood(vnHood3d, x, y, z);
                            int ct_max = MapHood(vnHood3d, x, y, z);
                            if (ct == ct_max) {
                                oxygen.Set(x, y, z, initial_oxygen);
                            }
                        }
                    }
                }
            }
                oxygen.Update();
                oxygen.DiffusionADI(D_O2);
            oxygen.Update();
        }
        }


    // FUNCTION TO DIFFUSE CARBON SUBSTRATES
    public void Carbon_Diffuse(double initial_glucose,double D_G) {
        for (int i = 0; i < 5; i++) {
            for (int x = 0; x < xDim; x++) {
                for (int z = 0; z < zDim; z++) {
                    for (int y = 0; y < yDim; y++) {
                        if (GetAgent(x, y, z) == null) {
                            int[] vnHood3d = Util.MooreHood3D(false);
                            int ct = MapEmptyHood(vnHood3d, x, y, z);
                            int ct_max = MapHood(vnHood3d, x, y, z);
                            if (ct == ct_max) {
                                glucose.Set(x, y, z, initial_glucose);
                            }
                        }
                          }
                    }
                }
                glucose.Update();
                glucose.DiffusionADI(D_G);
            }
            glucose.Update();
        }

    // FUNCTION TO PERFORM AGENT FUNCTIONS
    public void StepCells() {
        CleanAgents(); // CLEAN DEAD AGENTS IF NECESSARY
        // ITERATE OVER EACH AGENT
        for (Cell3D cell : this) {
            if (cell.metabolic_state != 10) { // IF THE CELL IS NOT DEAD
                int birth_counter = cell.Biomass_Divide(max_cell_index); // AGENT GROWTH AND DIVISION
                //  RESET MAXIMUM CELL INDEX (USED TO ASSIGN NEW CELL INDEX TO DAUGHTER CELLS) IF DIVISION EVENT OCCURED
                if (birth_counter > 0) {
                    max_cell_index = max_cell_index + 1;
                }

                // PERFORM CELL SHOVING
                for (int i = 0; i < 5; i++) {
                    cell.CalcMove();
                    cell.MoveDiv();
                }

            }
        }
            CleanAgents(); // CLEAN DEAD AGENTS IF NECESSARY
    }

    // FUNCTION TO RANDOMLY INITIALIZE AGENTS
    public void Initialize_Random(int b0) {
        max_cell_index = b0; // SET MAXIMUM CELL INDEX
        for (int j = 0; j <= max_cell_index-1; j++) { // j DEFINES THE AGENT INDEX
            double x = rn.Double(xDim); // X-COORDINATE
            double y = rn.Double(yDim); // Y-COORDINATE
            Random r = new Random();
            double old_mass =  1e-12 + (2e-12 - 1e-12) * r.nextDouble(); // SET AGENT BIOMASS
            int angle = rn.Int(360); // INITIALIZE CELL ANGLE
            int metabolic_state = 0;
            NewAgentPT(x, y, 0).Init(old_mass, angle, metabolic_state, j, 0, 0);
        }
    }

    // FUNCTION TO CALL MAIN METHOD
    static int count = 0;
    public static void mainCaller() throws IOException, NoSuchFieldException, IllegalAccessException {
        count++;
        while (count == 1) {
            main(null);
        }
    }

    public static void main(String[] args) throws IOException, NoSuchFieldException, IllegalAccessException {
        mainCaller(); // MAIN METHOD USED EXTERNALLY IN JAVA GATEWAY FILE
    }
}