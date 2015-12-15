#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <nlopt.h>

// uses core separations; modifications made to which minimized runs are excluded
// additional data generated from previous params in lammps for iterative fitting, then PARTIALLY OPTIMIZED in Gaussian (accesses data from HSEHsol_xyz2)

#define K_COUL 332.06342 // in kcal/mol, Angstroms, elementary charges
#define N_DOT_SIZES 6
#define N_RUNS 1450 // total number of runs across all dot sizes
#define N_ATOM_TYPES 9
#define N_ATOMS_TOTAL 179010 //all xyz points in all xyz files
// #define F_TOL 3.557 // tolerance on the force on one complex in kcal/mol/A based on Gaussian convergence criteria for Opt

#define PbS 0
#define S   1
#define PbO 2
#define OO	3
#define OH  4
#define HO  5
#define COO 6
#define CCO 7
#define HC  8

#define N_PARAMS 104

const int N_RUNS_BY_SIZE[] = {223, 292, 344, 323, 134, 134};
const double MASS_BY_TYPE[] = {207.2, 32.065, 207.2, 15.9994, 15.9994, 1.0079, 12.0107, 12.0107, 1.0079};


//for storing an energy and a force
typedef struct {
	double e, fx, fy, fz;
} t_E_and_F;

//calculate Morse, LJ, and Coulombic energy and force
t_E_and_F morse_lj_coul2(double D, double alpha, double r0, double epsilon, double sigma, double q1, double q2, double rx, double ry, double rz) {
    double r = sqrt(rx*rx + ry*ry + rz*rz);
   
    double morse_energy = 0;
    double morse_force = 0;
    double lj_energy = 0;
    double lj_force = 0;
   
    if(D!=0.0 && r<10.0) {
        double ex = exp(-alpha*(r-r0));
        morse_energy = D*( ex*ex - 2*ex );
        morse_force = 2*alpha*D*( ex*ex - ex );
    }
    
    if(epsilon!=0.0 && r<10.0) {
        double s6 = pow(sigma, 6);
        double r6 = pow(r, 6);
        lj_energy = 4*epsilon*( (s6/r6)*(s6/r6) - (s6/r6) );
        lj_force = -24*epsilon*s6*( r6 - 2*s6 )/(r6*r6*r);
    }
	
    double coul_energy = (K_COUL * q1 * q2)/r;
    double coul_force = (K_COUL * q1 * q2)/(r*r);
   
    t_E_and_F result;
    result.e = morse_energy + lj_energy + coul_energy;
    result.fx = (morse_force + lj_force + coul_force) * rx/r;
    result.fy = (morse_force + lj_force + coul_force) * ry/r;
    result.fz = (morse_force + lj_force + coul_force) * rz/r;
   
    return result;
}

//store an atom
typedef struct {
	int type;
	double x,y,z;
} t_atom;

//store a run (multiple atoms, energy)
typedef struct {
	int n_atoms_in_complex, n_atoms_in_dot;
	t_atom *atoms_in_complex;
	t_atom *atoms_in_dot;
	char optimized;
} t_run;

t_run runs[N_DOT_SIZES][344]; //all run data

//global params
double morse_D     [N_ATOM_TYPES][N_ATOM_TYPES];
double morse_alpha [N_ATOM_TYPES][N_ATOM_TYPES];
double morse_r0    [N_ATOM_TYPES][N_ATOM_TYPES];
double lj_epsilon  [N_ATOM_TYPES][N_ATOM_TYPES];
double lj_sigma    [N_ATOM_TYPES][N_ATOM_TYPES];
double coul_charge [N_ATOM_TYPES];

//set the global parameters according to the current values
void set_morse_lj_coul(const double *params) {

// Pb-dot - OH oxygen
	morse_D     [PbS][OH] = params[0];
	morse_alpha [PbS][OH] = params[1];
	morse_r0    [PbS][OH] = params[2];
// Pb-dot - Pb-complex
	morse_D     [PbS][PbO] = params[3];
	morse_alpha [PbS][PbO] = params[4];
	morse_r0    [PbS][PbO] = params[5];
// S - OH oxygen
	morse_D     [S][OH] = params[6];
	morse_alpha [S][OH] = params[7];
	morse_r0    [S][OH] = params[8];
// S - Pb-complex
	morse_D     [S][PbO] = params[9];
	morse_alpha [S][PbO] = params[10];
	morse_r0    [S][PbO] = params[11];
// Pb-dot - COO oxygen
	morse_D     [PbO][OO] = params[12];
	morse_alpha [PbO][OO] = params[13];
	morse_r0    [PbO][OO] = params[14];
// S - COO oxygen
	morse_D     [S][OO] = params[15];
	morse_alpha [S][OO] = params[16];
	morse_r0    [S][OO] = params[17];
// Pb-complex - OH oxygen
	morse_D     [PbO][OH] = params[18];
	morse_alpha [PbO][OH] = params[19];
	morse_r0    [PbO][OH] = params[20];
// Pb-complex - Pb-complex
	morse_D     [PbO][PbO] = params[21];
	morse_alpha [PbO][PbO] = params[22];
	morse_r0    [PbO][PbO] = params[23];
// OH oxygen - OH oxygen
	morse_D     [OH][OH] = params[24];
	morse_alpha [OH][OH] = params[25];
	morse_r0    [OH][OH] = params[26];
// Pb-complex - COO oxygen
	morse_D     [PbO][OO] = params[27];
	morse_alpha [PbO][OO] = params[28];
	morse_r0    [PbO][OO] = params[29];
// COO oxygen - OH oxygen
	morse_D     [OO][OH] = params[30];
	morse_alpha [OO][OH] = params[31];
	morse_r0    [OO][OH] = params[32];
// COO oxygen - COO oxygen
	morse_D     [OO][OO] = params[33];
	morse_alpha [OO][OO] = params[34];
	morse_r0    [OO][OO] = params[35];
// COO oxygen - COO carbon (+)
	morse_D     [OO][COO] = params[36];
	morse_alpha [OO][COO] = params[37];
	morse_r0    [OO][COO] = params[38];
// COO oxygen - CCO carbon (-)
	morse_D     [OO][CCO] = params[39];
	morse_alpha [OO][CCO] = params[40];
	morse_r0    [OO][CCO] = params[41];
// OH oxygen - COO carbon (+)
	morse_D     [OH][COO] = params[42];
	morse_alpha [OH][COO] = params[43];
	morse_r0    [OH][COO] = params[44];
// OH oxygen - CCO carbon (-)
	morse_D     [OH][CCO] = params[45];
	morse_alpha [OH][CCO] = params[46];
	morse_r0    [OH][CCO] = params[47];
// PbS - COO (+)
	morse_D     [PbS][COO] = params[48];
	morse_alpha [PbS][COO] = params[49];
	morse_r0    [PbS][COO] = params[50];
// PbS - HC
	//lj_epsilon [PbS][HC] = params[37];
	//lj_sigma   [PbS][HC] = params[38];
// PbO - COO (+)
	morse_D     [PbO][COO] = params[51];
	morse_alpha [PbO][COO] = params[52];
	morse_r0    [PbO][COO] = params[53];
// PbO - HC
	//lj_epsilon [PbO][HC] = params[45];
	//lj_sigma   [PbO][HC] = params[46];
// S - COO (+)
	morse_D     [S][COO] = params[54];
	morse_alpha [S][COO] = params[55];
	morse_r0    [S][COO] = params[56];
// PbS - CCO (-)
	morse_D     [PbS][CCO] = params[57];
	morse_alpha [PbS][CCO] = params[58];
	morse_r0    [PbS][CCO] = params[59];
// PbO - CCO (-)
	morse_D     [PbO][CCO] = params[60];
	morse_alpha [PbO][CCO] = params[61];
	morse_r0    [PbO][CCO] = params[62];
// S - CCO (-)
	morse_D     [S][CCO] = params[63];
	morse_alpha [S][CCO] = params[64];
	morse_r0    [S][CCO] = params[65];
// COO carbon - COO carbon
	morse_D     [COO][COO] = params[66];
	morse_alpha [COO][COO] = params[67];
	morse_r0    [COO][COO] = params[68];
// PbS - HO
	morse_D     [PbS][HO] = params[69];
	morse_alpha [PbS][HO] = params[70];
	morse_r0    [PbS][HO] = params[71];
// PbO - HO
	morse_D     [PbO][HO] = params[72];
	morse_alpha [PbO][HO] = params[73];
	morse_r0    [PbO][HO] = params[74];
// S - HO
	morse_D     [S][HO] = params[75];
	morse_alpha [S][HO] = params[76];
	morse_r0    [S][HO] = params[77];
// OH - HO
	morse_D     [OH][HO] = params[78];
	morse_alpha [OH][HO] = params[79];
	morse_r0    [OH][HO] = params[80];
// OO - HO
	morse_D     [OO][HO] = params[81];
	morse_alpha [OO][HO] = params[82];
	morse_r0    [OO][HO] = params[83];
// COO - HO
	morse_D     [COO][HO] = params[84];
	morse_alpha [COO][HO] = params[85];
	morse_r0    [COO][HO] = params[86];
// COO - CCO
	morse_D     [COO][CCO] = params[87];
	morse_alpha [COO][CCO] = params[88];
	morse_r0    [COO][CCO] = params[89];
//PbS - PbS
	morse_D     [PbS][PbS] = params[90];
	morse_alpha [PbS][PbS] = params[91];
	morse_r0    [PbS][PbS] = params[92];
// S - S
	morse_D     [S][S] = params[93];
	morse_alpha [S][S] = params[94];
	morse_r0    [S][S] = params[95];
// PbS - S
	morse_D     [PbS][S] = params[96];
	morse_alpha [PbS][S] = params[97];
	morse_r0    [PbS][S] = params[98];
// S - HC
	//lj_epsilon [S][HC] = params[53];
	//lj_sigma   [S][HC] = params[54];
	
	coul_charge[PbS] = params[99];
	coul_charge[S] = -coul_charge[PbS];
	coul_charge[PbO] = params[100];
	coul_charge[OH] = params[101];
	coul_charge[OO] = params[102];
	coul_charge[HO] = params[103];
	coul_charge[COO] = -(coul_charge[PbO] + coul_charge[OH] + coul_charge[OO] + coul_charge[OO] + coul_charge[HO] + coul_charge[CCO] + coul_charge[HC] + coul_charge[HC] + coul_charge[HC]);
	
	//Set params for t1>t2. assumes all above params [t1][t2] have been set such that t1<t2
	int t1,t2;
	for(t2=0; t2<N_ATOM_TYPES; t2++) {
		for(t1=t2+1; t1<N_ATOM_TYPES; t1++) {
			morse_D     [t1][t2] = morse_D     [t2][t1];
			morse_alpha [t1][t2] = morse_alpha [t2][t1];
			morse_r0    [t1][t2] = morse_r0    [t2][t1];
			lj_epsilon	[t1][t2] = lj_epsilon  [t2][t1];
			lj_sigma	[t1][t2] = lj_sigma    [t2][t1];
		}
	}
}

double *target_values;
int count_function_calls = 0;
FILE *parameter_log;
clock_t starting_clock;
double best_parameters[N_PARAMS];
double best_energy_error = 1e20;
double best_force_error = 1e20;
double best_torque_error = 1e20;
double best_energy_sum = 1e20;
double best_force_sum = 1e20;

//the function to be optimized
double error_function(unsigned params_count, const double *params, double *gradient, void *optional_data) {
	set_morse_lj_coul(params);
	double force_error = 0.0;
	double force_error_sum = 0.0;
	double energy_error = 0.0;
	double max_error = 1e-20;
	int max_index = 0;
	double energy_error_sum = 0.0;
	double torque_error = 0.0;
	int run_index=0, i, j, dot_size, run_number;
	for(dot_size=0; dot_size<N_DOT_SIZES; dot_size++) {
		for(run_number=0; run_number<N_RUNS_BY_SIZE[dot_size]; run_number++) {
			t_run run = runs[dot_size][run_number];
			int n_atoms_in_complex = run.n_atoms_in_complex, n_atoms_in_dot = run.n_atoms_in_dot;
			double comx = 0.0;
			double comy = 0.0;
			double comz = 0.0;
			double m_tot = 0.0;
			// calculate center of mass for one complex
			for(i=0; i<n_atoms_in_complex; i++) {
				t_atom c = run.atoms_in_complex[i];
				comx  += c.x*MASS_BY_TYPE[ c.type ];
				comy  += c.y*MASS_BY_TYPE[ c.type ];
				comz  += c.z*MASS_BY_TYPE[ c.type ];
				m_tot += 	 MASS_BY_TYPE[ c.type ];
			}
			comx = comx/m_tot;
			comy = comy/m_tot;
			comz = comz/m_tot;
			double ee = 0.0;
			double fx = 0.0;
			double fy = 0.0;
			double fz = 0.0;
			double tx = 0.0;
			double ty = 0.0;
			double tz = 0.0;
			for(i=0; i<n_atoms_in_complex; i++) {
				t_atom a = run.atoms_in_complex[i];
				for(j=0; j<n_atoms_in_dot; j++) {
					t_atom b = run.atoms_in_dot[j];
				
					double rx = b.x-a.x;
					double ry = b.y-a.y;
					double rz = b.z-a.z;
				
					double D =     morse_D[ a.type ][ b.type ];
					double alpha = morse_alpha[ a.type ][ b.type ];
					double r0 =    morse_r0[ a.type ][ b.type ];
				
					double epsilon = lj_epsilon[ a.type ][ b.type ];
					double sigma =   lj_sigma[ a.type ][ b.type ];
				
					double q1 = coul_charge[ a.type ];
					double q2 = coul_charge[ b.type ];
				
					t_E_and_F new_ef = morse_lj_coul2(D, alpha, r0, epsilon, sigma, q1, q2, rx, ry, rz);
					ee += new_ef.e;
					fx += new_ef.fx;
					fy += new_ef.fy;
					fz += new_ef.fz;
				}
				// calculate torque on one atom of the complex
				tx += fx*(a.x-comx);
				ty += fy*(a.y-comy);
				tz += fz*(a.z-comz);
			}

			if((run_index != 60) && (run_index != 65) && (run_index != 283) && (run_index != 285) && (run_index != 287) && (run_index != 290) && (run_index != 922) && (run_index != 925) && (run_index != 929) && (run_index != 930) && (run_index != 76) && (run_index != 85) && (run_index != 94) && (run_index != 103) && (run_index != 112) && (run_index != 121) && (run_index != 124) && (run_index != 130) && (run_index != 139) && (run_index != 142) && (run_index != 148) && (run_index != 157) && (run_index != 163) && (run_index != 169) && (run_index != 193) && (run_index != 412) && (run_index != 423) && (run_index != 456) && (run_index != 782) && (run_index != 786) && (run_index != 795) && (run_index != 801) && (run_index != 197) && (run_index != 199) && (run_index != 463) && (run_index != 464) && (run_index != 476) && (run_index != 487) && (run_index != 488) && (run_index != 489) && (run_index != 499) && (run_index != 510) && (run_index != 822) && (run_index != 823) && (run_index != 827) && (run_index != 832) && (run_index != 834) && (run_index != 839) && (run_index != 840) && (run_index != 853) && (run_index != 857) && (run_index != 1113) && (run_index != 1120) && (run_index != 1169) && (run_index != 1180) && (run_index != 168) && (run_index != 172) && (run_index != 194) && (run_index != 414) && (run_index != 417) && (run_index != 432) && (run_index != 437) && (run_index != 439) && (run_index != 449) && (run_index != 459) && (run_index != 462) && (run_index != 471) && (run_index != 766) && (run_index != 776) && (run_index != 781) && (run_index != 804) && (run_index != 1121) && (run_index != 1155) && ((run_index > 2 && run_index < 223) || (run_index > 225 && run_index < 291) || (run_index > 293 && run_index < 352) || (run_index > 410 && run_index < 515) || (run_index > 517 && run_index < 585) || (run_index > 587 && run_index < 646) || (run_index > 704 && run_index < 706) || (run_index > 764 && run_index < 859) || (run_index > 861 && run_index < 931) || (run_index > 933 && run_index < 992) || (run_index > 1050 && run_index < 1052) || (run_index > 1110 && run_index < 1182) || (run_index > 1184 && run_index < 1256) || (run_index > 1258 && run_index < 1316) || (run_index > 1318 && run_index < 1376) || (run_index > 1378))) { //total = 64+21+5-5+3+4+2-5+21-1+18-5=122 tentative - 16
				double e_dif = target_values[run_index] - ee;
				energy_error += e_dif*e_dif;
				double abs_e_dif = fabs(e_dif);
				energy_error_sum += abs_e_dif;
				if(abs_e_dif > max_error) {
					max_error = abs_e_dif;
					max_index = run_index;
				}
			
				if(run.optimized) { //only add force to results if run was optimized
					double f_dif = fx*fx + fy*fy + fz*fz;
					force_error  += f_dif * 0.3;
					force_error_sum += sqrt(f_dif);
					torque_error += (tx*tx + ty*ty + tz*tz) * 1e-4;
				}
			}
			run_index++;
		}
	}
	count_function_calls += 1;
	if(energy_error+force_error+torque_error < best_energy_error+best_force_error+best_torque_error) { //store best parameters
		memcpy(best_parameters, params, N_PARAMS*sizeof(double));
		best_energy_error = energy_error;
		best_force_error = force_error;
		best_torque_error = torque_error;
		best_energy_sum = energy_error_sum;
		best_force_sum = force_error_sum;
	}
	if( count_function_calls%1000 == 0 ) {
        for(i=0; i<N_PARAMS; ++i) {
            fprintf(parameter_log, "%.7g ", best_parameters[i]);
        }
        printf("Best energy error: %.7g\nBest force error: %.7g\nBest torque error: %.7g\nEnergy Error Sum: %.7g\nForce Error Sum: %.7g\nMax Error: %.7g at index %d\n", best_energy_error, best_force_error, best_torque_error, best_energy_sum, best_force_sum, max_error, max_index);
        fprintf(parameter_log, "\nEnergy Error = %g, Force Error = %g, Torque Error = %g, Time = %gs\n", best_energy_error, best_force_error, best_torque_error, ((double)clock()-starting_clock)/CLOCKS_PER_SEC );
        fflush(parameter_log);
    }
	return force_error + energy_error + torque_error;
}

int main( int argc, char *argv[] ) {
	t_atom *atoms = (t_atom*) malloc( N_ATOMS_TOTAL * sizeof(t_atom) ); //Put each new atom here. Put references to this in run.atom for each run. 
	
	printf("Allocating %ld bytes\n", N_ATOMS_TOTAL * sizeof(t_atom));
	
	target_values = (double*) malloc( (N_RUNS) * sizeof(double) ); //want to match N_RUNS energies
	
	int target_value_index = 0;
	int atom_index = 0;
	int sum_atom_counts = 0;
	char filename[100];
	int dot_size, i;
	int fscanf_result;
	for(dot_size=0; dot_size<N_DOT_SIZES; dot_size++) {
		for(i=0; i<N_RUNS_BY_SIZE[dot_size]; i++) { //read through all values of i
			//atoms in dot
			sprintf(filename, "/fs/home/ja539/build/levmar-2.6/HSEHsol_xyz2/HSEHsol_dot%d_%d.xyz", dot_size+1, i);
			FILE *f = fopen(filename, "r");
			if( f == NULL ) { break; } //no larger values of i for this dot size
			
			fscanf(f, "%d", &runs[dot_size][i].n_atoms_in_dot); //read n atoms
			fscanf(f, "%lf %hhd", &target_values[target_value_index], &runs[dot_size][i].optimized); //read energy and whether the run is optimized
			
			runs[dot_size][i].atoms_in_dot = &atoms[atom_index];
			sum_atom_counts += runs[dot_size][i].n_atoms_in_dot; //debugging

			int type;
			double x, y, z;
			while( (fscanf_result=fscanf(f, "%d %lf %lf %lf", &type, &x, &y, &z)) > 0 ) {
				atoms[atom_index].type = type;
				atoms[atom_index].x = x;
				atoms[atom_index].y = y;
				atoms[atom_index].z = z;
				if( atom_index >= sum_atom_counts ) {
					printf("Bad atom %d/%d: %d %f %f %f\n", atom_index, sum_atom_counts, type, x, y, z);
					exit(1);
				}
				atom_index++;
			}
			fclose(f);
			
			//atoms in complex
			sprintf(filename, "/fs/home/ja539/build/levmar-2.6/HSEHsol_xyz2/HSEHsol_complex%d_%d.xyz", dot_size+1, i);
			//puts(filename); //debugging
			f = fopen(filename, "r");
			if( f == NULL ) { puts("Error: this dot has no complex."); exit(1); }
			
			fscanf(f, "%d", &runs[dot_size][i].n_atoms_in_complex); //read n atoms
			double eng; char opt;
			fscanf(f, "%lf %hhd", &eng, &opt);

			if(eng!=target_values[target_value_index] || opt!=runs[dot_size][i].optimized) {
				printf("Inconsistent complex and dot: %f %f, %d %d\n", eng, target_values[target_value_index], opt, runs[dot_size][i].optimized);
				exit(1);
			}
			
			runs[dot_size][i].atoms_in_complex = &atoms[atom_index];
			sum_atom_counts += runs[dot_size][i].n_atoms_in_complex; //debugging;
			target_value_index++;
			
			while( fscanf(f, "%d %lf %lf %lf", &type, &x, &y, &z) > 0 ) {
				atoms[atom_index].type = type;
				atoms[atom_index].x = x;
				atoms[atom_index].y = y;
				atoms[atom_index].z = z;
				if( atom_index >= sum_atom_counts ) {
					printf("Bad atom %d/%d: %d %f %f %f\n", atom_index, sum_atom_counts, type, x, y, z);
					exit(1);
				}
				atom_index++;
			}
			fclose(f);
			
		}
	}
	
	puts("Done reading atoms");
	
	//set up global parameters
	double epsilon[N_ATOM_TYPES]; //use to store i-i interactions before setting i-j
	for(i=0; i<N_ATOM_TYPES; i++) {
		epsilon[i]=0.0;
	}
	double sigma[N_ATOM_TYPES];
	
	//set according to OPLS
	//coul_charge[OO] = 0.0; treated as variable
	//coul_charge[OH] = 0.0; treated as variable
	//coul_charge[HO] = 0.3; treated as variable
	//coul_charge[COO] = 0.7;treated as variable
	coul_charge[CCO] = -0.28;
	coul_charge[HC] = 0.06;
	
	epsilon[OO] = 0.21;
	sigma  [OO] = 2.96;
	epsilon[OH] = 0.25;
	sigma  [OH] = 3.2;
	epsilon[HO] = 0.0;
	sigma  [HO] = 0.0;
	epsilon[COO] = 0.1050;
	sigma  [COO] = 3.75;
	epsilon[CCO] = 0.0660;
	sigma  [CCO] = 3.5;
	epsilon[HC] = 0.03;
	sigma  [HC] = 2.5;
	
	int t1,t2;
	for(t1=0; t1<N_ATOM_TYPES; t1++) {
		for(t2=0; t2<N_ATOM_TYPES; t2++) {
			morse_D     [t1][t2] = 0.0; //zero energy
			morse_alpha [t1][t2] = 1.0; //dummy param
			morse_r0    [t1][t2] = 3.0; //dummy param
			
			if( epsilon[t1]==0.0 || epsilon[t2]==0.0 ) {
				lj_epsilon[t1][t2] = 0.0; //zero energy
				lj_sigma  [t1][t2] = 3.0; //dummy param
			}
			else {
				lj_epsilon[t1][t2] = sqrt( epsilon[t1]*epsilon[t2] ); //OPLS mixing rules
				lj_sigma  [t1][t2] = sqrt( sigma[t1]*sigma[t2] ); //OPLS mixing rules
			}
		}
	}
	
	lj_epsilon[OO][OO] = 0.0;
	lj_epsilon[COO][COO] = 0.0;
	lj_epsilon[COO][CCO] = 0.0;
	lj_epsilon[CCO][COO] = 0.0;
	lj_epsilon[OH][OO] = 0.0;
	lj_epsilon[OO][OH] = 0.0;
	lj_epsilon[OH][OH] = 0.0;
	lj_epsilon[COO][OH] = 0.0;
	lj_epsilon[OH][COO] = 0.0;
	lj_epsilon[COO][OO] = 0.0;
	lj_epsilon[OO][COO] = 0.0;
	lj_epsilon[CCO][OO] = 0.0;
	lj_epsilon[OO][CCO] = 0.0;
	lj_epsilon[OH][CCO] = 0.0;
	lj_epsilon[CCO][OH] = 0.0;
	lj_epsilon[COO][CCO] = 0.0;
	lj_epsilon[CCO][COO] = 0.0;
	
	//additional changes to increase stability:
	//additional changes to increase stability:
	//lj_epsilon[OH][CCO]  = 0.01;
	//lj_epsilon[CCO][OH]  = 0.01;
	lj_epsilon[OH][HC]   = 0.01;
	lj_epsilon[HC][OH]   = 0.01;
	lj_epsilon[OO][HC]   = 0.01;
	lj_epsilon[HC][OO]   = 0.01;
	//lj_epsilon[OO][HO]   = 0.01;
	//lj_epsilon[HO][OO]   = 0.01;
	//lj_epsilon[OO][CCO]  = 0.01;
	//lj_epsilon[CCO][OO]  = 0.01;
	//lj_epsilon[PbS][COO] = 0.01;
	//lj_epsilon[PbS][CCO] = 0.01;
	//lj_epsilon[PbS][HO]  = 0.01;
	lj_epsilon[PbS][HC]  = 0.01;
	//lj_epsilon[COO][PbS] = 0.01;
	//lj_epsilon[CCO][PbS] = 0.01;
	//lj_epsilon[HO][PbS]  = 0.01;
	lj_epsilon[HC][PbS]  = 0.01;
	//lj_epsilon[PbO][COO] = 0.01;
	//lj_epsilon[PbO][CCO] = 0.01;
	//lj_epsilon[PbO][HO]  = 0.01;
	lj_epsilon[PbO][HC]  = 0.01;
	//lj_epsilon[COO][PbO] = 0.01;
	//lj_epsilon[CCO][PbO] = 0.01;
	//lj_epsilon[HO][PbO]  = 0.01;
	lj_epsilon[HC][PbO]  = 0.01;
	//lj_epsilon[S][COO]   = 0.01;
	//lj_epsilon[S][CCO]   = 0.01;
	//lj_epsilon[S][HO]    = 0.01;
	lj_epsilon[S][HC]    = 0.01;
	//lj_epsilon[COO][S]   = 0.01;
	//lj_epsilon[CCO][S]   = 0.01;
	//lj_epsilon[HO][S]    = 0.01;
	lj_epsilon[HC][S]    = 0.01;
	
	//lj_epsilon [S][OH] = 0.01;
	//lj_epsilon [OH][S] = 0.01;
	//lj_sigma   [S][OH] = 3.0;
	//lj_sigma   [OH][S] = 3.0;
	//lj_epsilon [S][OO] = 0.01;
	//lj_epsilon [OO][S] = 0.01;
	//lj_sigma   [S][OO] = 3.0;
	//lj_sigma   [OO][S] = 3.0;
	//lj_epsilon [OO][OH] = 0.01;
	//lj_epsilon [OH][OO] = 0.01;
	//lj_sigma   [OO][OH] = 3.0;
	//lj_sigma   [OH][OO] = 3.0;
	
	//lj_sigma[OH][CCO]  = 4.0;
	//lj_sigma[CCO][OH]  = 4.0;
	lj_sigma[OH][HC]   = 2.5;
	lj_sigma[HC][OH]   = 2.5;
	lj_sigma[OO][HC]   = 2.5;
	lj_sigma[HC][OO]   = 2.5;
	//lj_sigma[OO][HO]   = 2.5;
	//lj_sigma[HO][OO]   = 2.5;
	//lj_sigma[OO][CCO]  = 4.0;
	//lj_sigma[CCO][OO]  = 4.0;
	//lj_sigma[PbS][COO] = 4.0;
	//lj_sigma[PbS][CCO] = 4.0;
	//lj_sigma[PbS][HO]  = 2.5;
	lj_sigma[PbS][HC]  = 2.5;
	//lj_sigma[COO][PbS] = 4.0;
	//lj_sigma[CCO][PbS] = 4.0;
	//lj_sigma[HO][PbS]  = 2.5;
	lj_sigma[HC][PbS]  = 2.5;
	//lj_sigma[PbO][COO] = 4.0;
	//lj_sigma[PbO][CCO] = 4.0;
	//lj_sigma[PbO][HO]  = 2.5;
	lj_sigma[PbO][HC]  = 2.5;
	//lj_sigma[COO][PbO] = 4.0;
	//lj_sigma[CCO][PbO] = 4.0;
	//lj_sigma[HO][PbO]  = 2.5;
	lj_sigma[HC][PbO]  = 2.5;
	//lj_sigma[S][COO]   = 4.0;
	//lj_sigma[S][CCO]   = 4.0;
	//lj_sigma[S][HO]    = 2.5;
	lj_sigma[S][HC]    = 2.5;
	//lj_sigma[COO][S]   = 4.0;
	//lj_sigma[CCO][S]   = 4.0;
	//lj_sigma[HO][S]    = 2.5;
	lj_sigma[HC][S]    = 2.5;
	//read mutable parameters from text file
	
	double *parameters = (double*) malloc( N_PARAMS * sizeof(t_run) );
	double *lower_bounds = (double*) malloc( N_PARAMS * sizeof(t_run) );
	double *upper_bounds = (double*) malloc( N_PARAMS * sizeof(t_run) );
	
	FILE *f = fopen(argv[1], "r");
	
	double param, lower, upper;
	
	i = 0;
	while( fscanf(f, "%lf %lf %lf", &param, &lower, &upper) != EOF ) {
		parameters[i] = param;
		lower_bounds[i] = lower;
		upper_bounds[i] = upper;
		i++;
	}
	
	set_morse_lj_coul(parameters);
	parameter_log = fopen(argv[4],"a");
	puts("Running optimization...");
	
	nlopt_opt opt;
	if(argv[5][0] == '1') { //MLSL-LDS  (deterministic)
		opt = nlopt_create(NLOPT_G_MLSL_LDS, N_PARAMS);
		nlopt_opt local_opt = nlopt_create(NLOPT_LN_SBPLX, N_PARAMS);
		nlopt_set_local_optimizer(opt, local_opt);
	} else if(argv[5][0] == '2') { //MLSL
		opt = nlopt_create(NLOPT_G_MLSL, N_PARAMS);
		nlopt_opt local_opt = nlopt_create(NLOPT_LN_SBPLX, N_PARAMS);
		nlopt_set_local_optimizer(opt, local_opt);
	} else if(argv[5][0] == '3') { //DIRECT  (deterministic)
		opt = nlopt_create(NLOPT_GN_DIRECT, N_PARAMS);
	} else if(argv[5][0] == '4') { //DIRECT-L  (deterministic)
		opt = nlopt_create(NLOPT_GN_DIRECT_L, N_PARAMS);
	} else if(argv[5][0] == '5') { //CRS
		opt = nlopt_create(NLOPT_GN_CRS2_LM, N_PARAMS);
	} else if(argv[5][0] == '6') { //ISRES  (deterministic)
		opt = nlopt_create(NLOPT_GN_ISRES, N_PARAMS);
	} else if(argv[5][0] == '7') { //COBYLA
		opt = nlopt_create(NLOPT_LN_COBYLA, N_PARAMS);
	} else { // SBPLX (local)
		opt = nlopt_create(NLOPT_LN_SBPLX, N_PARAMS);
	}
	
	nlopt_set_lower_bounds(opt, lower_bounds);
	nlopt_set_upper_bounds(opt, upper_bounds);
	//int maxeval = 100000;
	//int maxeval = 1000;
	//int maxeval = 200000;
	int maxeval = 2000000;
	//int maxeval = 10000000;
	nlopt_set_maxeval(opt, maxeval);
	nlopt_set_min_objective(opt, error_function, NULL);
	
	double error; /* the minimum objective value, upon return */
	
	starting_clock = clock();
	
	int return_code = nlopt_optimize(opt, parameters, &error);
	if(return_code < 0) {
		printf("nlopt failed with code %d\n", return_code);
	}
	else {
		printf("Found minimum at\n");
		for(i=0; i<N_PARAMS; ++i) {
			printf("%.7g ", parameters[i]);
		}
		printf("\nError = %g.\n", error);
		printf("Error per run = %g.\n", error/N_RUNS);
		
		FILE *e = fopen(argv[2],"w");
		fprintf(e, "%f\n", error);
		fclose(e);

		// Write results to text file.
		FILE *h = fopen(argv[3],"w");
		for( i=0; i<N_PARAMS; i++) {
			fprintf(h, "%f %f %f\n", parameters[i], lower_bounds[i], upper_bounds[i]);
		}
	
		fclose(h);
	}
	
	return 0;
}
