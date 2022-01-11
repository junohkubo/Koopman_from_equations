#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//--- The following constants can be changed to obtain more accurate results while a longer computational time is needed.
#define MAX_STATE_RANGE 30 // Perform cutoff for the state over MAX_STATE_RANGE for each dimension.
#define MAX_M 50 // Perform the resolvent calculation up to MAX_M.
#define CONV_DICIMAL_PLACE 5 // Evaluate the result with convergence at the decimal place CONV_DICIMAL_PLACE.

#define MAX_ID 10000
#define MAX_LEN_STR 4096

#define TRUE 1
#define FALSE 0
typedef long long int int64; // 64 bit integer

// Declare global variables from .inp files.
// *** NOTE:
// *** ori_event: The original events in the input file
// *** event: Events with the same state change (but different rate constants)
static int NUM_ORI_EVENT; // The number of ori_events
static int NUM_EVENT; // The number of events (without state-change dupilication)
static int DIM_STATE; // The dimension of state (the number of components)
static int NUM_CONSTANT; // The number of constants in the input file
//--- The followings are used for the correspondence between ori_events and events.
static int **ORI_EVENT;
static int **EVENT;
static double *ORI_RATE;
static int **ORI_INDEX4RATE;
static int **ORI_INDEX4CONST;
static int *LIST_LEN_E2ORI;
static int *LIST_START_INDEX_E2ORI;
//--- The followings are needed to calculate entities of the resolvent.
static int ID_DIAGONAL_EVENT;

static double *PRUNING_COEFF_POS;
static double *PRUNING_COEFF_NEG;
static double FINAL_TIME;
static int *INITIAL_STATE;
static int *TARGET_STATE;
static double *CONSTS;


// Declare working global variables.
static double *ORI_COEFF;
static int **ORI_STATE;
static double *ORI_EVENT_CONST_RATE;
static double *COEFF;
static int **STATE;


// Declare subroutines for the evaluation.
void load_event_file(char *filename);
void load_parameter_file(char *filename);
void initialize_settings_for_dual_comp(void);
void set_pruning_coeffs(void);
void set_state_from_str(int *state, char *str);
//int is_pruned(int *current, int length);
int is_pruned(int *current, int *target, int length);
int is_out_of_range(int *current);
double evaluate_resolvent_for_1st(double *c, double dt, int *target);

double evaluate_resolvent(double dt);
void calculate_rate_coeffs(void);
double occur_event(int event, int *index, double dt_M);

void aitken(double *seq, int len);


// Declare subroutines for file open and memory allocations.
FILE *f_open(const char *path, const char *mode);
void *clear_alloc(size_t size, size_t nobj);
void **clear_alloc_2d( size_t size, size_t row, size_t column);
#define free_2d(p) {free(p[0]);free(p);}

void free_memories_allocated_in_load_event_file(void);


int main(int argc, char *argv[]){
  if(argc != 3){
    fprintf(stderr, "Usage: ./app_dual_computation event_input_file parameter_input_file\n" );
    fprintf(stderr, "example: ./analysis tmp_event.inp tmp_parameters.inp\n" );
    exit(1);
  }
  clock_t start_clock, end_clock;
  start_clock = clock();
  
  // Load the .inp files (event and parameter).
  load_event_file(argv[1]);
  INITIAL_STATE = (int*)clear_alloc(sizeof(int), DIM_STATE);
  TARGET_STATE = (int*)clear_alloc(sizeof(int), DIM_STATE);
  CONSTS = (double*)clear_alloc(sizeof(double), NUM_CONSTANT);
  load_parameter_file(argv[2]);

  // Initialize settings.
  initialize_settings_for_dual_comp();

  // Perform the main calculation.
  double comp = evaluate_resolvent(FINAL_TIME);
  printf("%e\n", comp);

  // Output the CPU time.
  end_clock = clock();
  double cpu_time;
  cpu_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
  fprintf(stderr, "CPU TIME: %f\n", cpu_time);

  // Free the allocated memories.
  free(ORI_EVENT_CONST_RATE);
  free_2d(ORI_STATE);
  free_2d(STATE);
  free(ORI_COEFF);
  free(COEFF);
  free_memories_allocated_in_load_event_file();
  free(PRUNING_COEFF_NEG);
  free(PRUNING_COEFF_POS);
  free(CONSTS);
  free(TARGET_STATE);
  free(INITIAL_STATE);
  
  return 0;  
}

// -----------------------------------------------------------
// Subroutines 
// -----------------------------------------------------------
void load_event_file(char *filename){
  FILE *fp;
  char buf[MAX_LEN_STR];
  char str_state_changes[MAX_LEN_STR];
  char str_state_changes_previous[MAX_LEN_STR];
  char str_rate_coeff[MAX_LEN_STR];
  char str_rate_factors[MAX_LEN_STR];
  char str_rate_constants[MAX_LEN_STR];
  fp = f_open(filename, "r");
  // -- Skip the first line.
  fgets(buf, sizeof(buf), fp);
  // -- Read parameters.
  fscanf(fp, "%d %d %d %d\n", &NUM_ORI_EVENT, &NUM_EVENT, &DIM_STATE, &NUM_CONSTANT);
  // -- Skip the third line.
  fgets(buf, sizeof(buf), fp);
  // Allocate memories.
  ORI_EVENT = (int**)clear_alloc_2d(sizeof(int), NUM_ORI_EVENT, DIM_STATE);
  EVENT = (int**)clear_alloc_2d(sizeof(int), NUM_EVENT, DIM_STATE);
  ORI_RATE = (double*)clear_alloc(sizeof(double), NUM_ORI_EVENT);
  ORI_INDEX4RATE = (int**)clear_alloc_2d(sizeof(int), NUM_ORI_EVENT, DIM_STATE);
  ORI_INDEX4CONST = (int**)clear_alloc_2d(sizeof(int), NUM_ORI_EVENT, NUM_CONSTANT);
  LIST_LEN_E2ORI = (int*)clear_alloc(sizeof(int), NUM_EVENT + 1); // '+1' is needed for the case in which the diagonal part is zero.
  LIST_START_INDEX_E2ORI = (int*)clear_alloc(sizeof(int), NUM_EVENT);
  for(int i = 0; i < NUM_EVENT+1; i++){
    LIST_LEN_E2ORI[i] = 0;
  }
  int e = 0;
  int ori_e = 0;
  int is_there_a_diagonal_event = FALSE;
  strcpy(str_state_changes_previous, "");
  while(fgets(buf, sizeof(buf), fp) != NULL){
    char str_tmp[4096];
    if(*buf == '\n' || buf[0] == ' ') continue;
    sscanf(buf,"%[^/] / %[^/] / %[^/] / %[^#]", str_state_changes, str_rate_coeff, str_rate_factors, str_rate_constants);
    {
      // Read state_change.
      if(strcmp(str_state_changes, str_state_changes_previous) != 0){
	// Set EVENT and ORI_EVENT.
	char tmp[4096];
	strcpy(tmp, str_state_changes);
	char *tp;
	int c = 0;
	int is_zero_change = TRUE;
	tp = strtok(tmp, " ,[]");
	ORI_EVENT[ori_e][c] = atoi(tp);
	EVENT[e][c] = atoi(tp);
	if(EVENT[e][c] != 0){
	  is_zero_change = FALSE;
	}
	c++;
	while(tp != NULL){
	  tp = strtok(NULL," ,[]");
	  if (tp != NULL){
	    ORI_EVENT[ori_e][c] = atoi(tp);
	    EVENT[e][c] = atoi(tp);
	    if(EVENT[e][c] != 0){
	      is_zero_change = FALSE;
	    }
	    c++;
	  }
	}
	// Set the correspondence from the event to the ori_event.
	LIST_START_INDEX_E2ORI[e] = ori_e;
	LIST_LEN_E2ORI[e]++;
	// Check whether there is an no-state-change event or not.
	if(is_zero_change == TRUE){
	  is_there_a_diagonal_event = TRUE;
	  ID_DIAGONAL_EVENT = e;
	}
	e++;
      }else{
	// Set only ORI_EVENT.
	char tmp[4096];
	strcpy(tmp, str_state_changes);
	char *tp;
	int c = 0;
	int is_zero_change = TRUE;
	tp = strtok(tmp, " ,[]");
	ORI_EVENT[ori_e][c] = atoi(tp);
	c++;
	while(tp != NULL){
	  tp = strtok(NULL," ,[]");
	  if (tp != NULL){
	    ORI_EVENT[ori_e][c] = atoi(tp);
	    c++;
	  }
	}
	// Set the correspondence from the event to the ori_event.
	LIST_LEN_E2ORI[e-1]++;
      }
      strcpy(str_state_changes_previous, str_state_changes);
      {
	// Read rate_coeff.
	ORI_RATE[ori_e] = atof(str_rate_coeff);
      }
      {
	// Read rate_factors.
	char *tp;
	int c = 0;
	tp = strtok(str_rate_factors, " ,[]");
	ORI_INDEX4RATE[ori_e][c] = atoi(tp);
	c++;
	while(tp != NULL){
	  tp = strtok(NULL," ,[]");
	  if (tp != NULL){
	    ORI_INDEX4RATE[ori_e][c] = atoi(tp);
	    c++;
	  }
	}
      }
      {
	// Read index for constant.
	char *tp;
	int c = 0;
	tp = strtok(str_rate_constants, " ,[]");
	ORI_INDEX4CONST[ori_e][c] = atoi(tp);
	c++;
	while(tp != NULL){
	  tp = strtok(NULL," ,[]");
	  if (tp != NULL){
	    ORI_INDEX4CONST[ori_e][c] = atoi(tp);
	    c++;
	  }
	}
      }
      ori_e++;
    }
  }
  fclose(fp);
  if(is_there_a_diagonal_event == FALSE){
    ID_DIAGONAL_EVENT = NUM_EVENT;
  }
}

void load_parameter_file(char *filename){
  FILE *fp;
  char buf[MAX_LEN_STR];
  fp = f_open(filename, "r");
  // -- Read the final time.
  fgets(buf, sizeof(buf), fp); // Skip a comment.
  fgets(buf, sizeof(buf), fp);
  sscanf(buf, "%lf ", &FINAL_TIME);
  // -- Read the initial state.
  fgets(buf, sizeof(buf), fp); // Skip a comment.
  fgets(buf, sizeof(buf), fp);
  set_state_from_str(INITIAL_STATE, buf);
  // -- Read the target state.
  fgets(buf, sizeof(buf), fp); // Skip a comment.
  fgets(buf, sizeof(buf), fp);
  set_state_from_str(TARGET_STATE, buf);
  // -- Read the parameters.
  fgets(buf, sizeof(buf), fp); // Skip a comment.
  for(int i=0; i<NUM_CONSTANT; i++){
    fgets(buf, sizeof(buf), fp);
    double comp;
    sscanf(buf, "%lf ", &comp);
    CONSTS[i] = comp;
  }
  fclose(fp);
}
  
void free_memories_allocated_in_load_event_file(void){
  free_2d(ORI_EVENT);
  free_2d(EVENT);
  free(ORI_RATE);
  free_2d(ORI_INDEX4RATE);
  free_2d(ORI_INDEX4CONST);
  free(LIST_LEN_E2ORI);
  free(LIST_START_INDEX_E2ORI);
}

void initialize_settings_for_dual_comp(void){
  /// Set coefficients for the pruning.
  set_pruning_coeffs();
  // Memory allocation.
  COEFF = (double*)clear_alloc(sizeof(double), MAX_ID);
  ORI_COEFF = (double*)clear_alloc(sizeof(double), MAX_ID);
  STATE = (int**)clear_alloc_2d(sizeof(int), MAX_ID, DIM_STATE);
  ORI_STATE = (int**)clear_alloc_2d(sizeof(int), MAX_ID, DIM_STATE);
  ORI_EVENT_CONST_RATE = (double*)clear_alloc(sizeof(double), NUM_ORI_EVENT);
}

void set_pruning_coeffs(void){
  // Find the maximum step in the positive and negative directions in each dimension.
  PRUNING_COEFF_POS = (double*)clear_alloc(sizeof(double), DIM_STATE);
  PRUNING_COEFF_NEG = (double*)clear_alloc(sizeof(double), DIM_STATE);
  for(int d = 0; d < DIM_STATE; d++){
    int max = 0;
    int min = 0;
    for(int e = 0; e < NUM_EVENT; e++){
      if(EVENT[e][d] > max){
	max = EVENT[e][d];
      }
      if(EVENT[e][d] < min){
	min = EVENT[e][d];
      }
    }
    if(max != 0){
      PRUNING_COEFF_POS[d] = 1.0/(double)max;
    }else{
      fprintf(stderr, "There is no event with positive direction in %d th dimension.", d);
    }
    if(min != 0){
      PRUNING_COEFF_NEG[d] = 1.0/(double)min;
    }else{
      fprintf(stderr, "There is no event with negative direction in %d th dimension.", d);
    }
  }
}

void set_state_from_str(int *state, char *str){
  char *tp;
  int d = 0;
  tp = strtok(str, " ,'\"");
  state[d] = atoi(tp);
  d++;
  while (tp != NULL){
    tp = strtok(NULL," ,'\"");
    if (tp != NULL){
      state[d] = atoi(tp);
      d++;
    }
  }
  if(d != DIM_STATE){
    fprintf(stderr, "The dimension of initial states is not adequate!\n");
    exit(1);
  }
}

int is_pruned(int *current, int *target, int length){
  // Calculate the number of steps required in each dimension and take the sum of them.
  double dist;
  for(int d = 0; d < DIM_STATE; d++){
    double comp = (double)(current[d]-target[d]);
    if(comp >= 0.0){
      double comp_pos = PRUNING_COEFF_POS[d] * comp;
      if(dist < comp_pos){
	dist = comp_pos;
      }
    }else{
      double comp_neg = PRUNING_COEFF_NEG[d] * comp;
      if(dist < comp_neg){
	dist = comp_neg;
      }
    }
  }
  // Judge whether the distance is over the length or not.
  if(dist <= length){
    return FALSE; // not to be pruned
  }else{
    return TRUE; // to be pruned
  }
}

int is_out_of_range(int *current){
  for(int d = 0; d < DIM_STATE; d++){
    if(current[d] > MAX_STATE_RANGE){
      return TRUE;
    }
  }
  return FALSE;
}

double evaluate_resolvent(double dt){
  int M;
  int64 num_of_terms, pre_num_of_terms;
  double pre_result;
  double current_result;
  double working_coeffs;
  int *working_state;
  working_state = (int*)clear_alloc(sizeof(int), DIM_STATE);
	    
  // Evalate constant rates (coefficients).
  for(int e = 0; e < NUM_ORI_EVENT; e++){
    double comp = ORI_RATE[e];
    for(int i = 0; i < NUM_CONSTANT; i++){
      comp = comp * pow(CONSTS[i], (double)ORI_INDEX4CONST[e][i]);
    }
    ORI_EVENT_CONST_RATE[e] = comp;
  }

  // Prepare sequences of differentials.
  // NOTE: To evaluate the extrapolation from the resolvent calculation, sequences of differentials are evaluated.
  //       The weighted sum is evaluated as the final result.
  int num_diff = 4; // Calculate the sequences of differentials up to the order num_diff.
  double **seq;
  seq = (double**)clear_alloc_2d(sizeof(double), num_diff, MAX_M+1);
  M = num_diff;
  
  // Prepare the temporal output file for the sequences of differentials.
  FILE *fp;
  fp = f_open("seq.dat", "w");

  int is_conv = FALSE;
  while(M <= MAX_M && is_conv == FALSE){
    double dt_M = dt / (double)M;
    double first_rate, final_rate;
    {
      // Calculate the first and final rates. These are multiplied in the final result.
      double rate = 0.0;
      for(int j = 0; j < LIST_LEN_E2ORI[ID_DIAGONAL_EVENT]; j++){
	int pre_e = LIST_START_INDEX_E2ORI[ID_DIAGONAL_EVENT] + j;
	double state_dependent_rate = 1.0;
	for(int d = 0; d < DIM_STATE; d++){
	  for(int fact = 0; fact < ORI_INDEX4RATE[pre_e][d]; fact++){
	    state_dependent_rate = state_dependent_rate * (double)(INITIAL_STATE[d] - fact);
	  }
	}
	rate = rate + state_dependent_rate * ORI_EVENT_CONST_RATE[pre_e];
      }
      first_rate = 1.0 - dt_M*rate;
      // Thr followings are not needed when the diagonal part does not include a constant.
      rate = 0.0;
      for(int j = 0; j < LIST_LEN_E2ORI[ID_DIAGONAL_EVENT]; j++){
	int pre_e = LIST_START_INDEX_E2ORI[ID_DIAGONAL_EVENT] + j;
	double state_dependent_rate = 1.0;
	for(int d = 0; d < DIM_STATE; d++){
	  for(int fact = 0; fact < ORI_INDEX4RATE[pre_e][d]; fact++){
	    state_dependent_rate = state_dependent_rate * (double)(TARGET_STATE[d] - fact);
	  }
	}
	rate = rate + state_dependent_rate * ORI_EVENT_CONST_RATE[pre_e];
      }
      final_rate = 1.0/(1.0 - dt_M*rate);
    }

    // Perform the initialization for each M.
    for(int d = 0; d < DIM_STATE; d++){
      STATE[0][d] = INITIAL_STATE[d];
    }
    COEFF[0] = 1.0;
    num_of_terms = 1;
    // Calculate coefficients and states iteratively.
    for(int depth = 1; depth <= M; depth++){
      // Copy the obtained states in the previous calculation.
      pre_num_of_terms = num_of_terms;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	ORI_COEFF[i] = COEFF[i];
	for(int d = 0; d < DIM_STATE; d++){
	  ORI_STATE[i][d] = STATE[i][d];
	}
      }
      num_of_terms = 0;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	for(int e = 0; e < NUM_EVENT; e++){
	  for(int d = 0; d < DIM_STATE; d++){
	    working_state[d] = ORI_STATE[i][d];
	  }
	  working_coeffs = occur_event(e, working_state, dt_M);
	  if(working_coeffs == 0.0){
	    continue;
	  }
	  // Check whether the new state is a candidate or not.
	  if(is_pruned(working_state, TARGET_STATE, M-depth+1) == TRUE){
	    continue;
	  }
	  if(is_out_of_range(working_state) == TRUE){
	    continue;
	  }
	  // Check whether the state has already appeared or not.
	  int is_already_appeared = FALSE;
	  int64 detected_index;
	  for(int64 j = 0; j < num_of_terms; j++){
	    int is_different = FALSE;
	    for(int d = 0; d < DIM_STATE; d++){
	      if(working_state[d] != STATE[j][d]){
		is_different = TRUE;
		break;
	      }
	    }
	    if(is_different == TRUE){
	      continue;
	    }else{
	      is_already_appeared = TRUE;
	      detected_index = j;
	      break;
	    }
	  }
	  if(is_already_appeared == TRUE){
	    // Revise the obtained coefficients.
	    // [CAUTION] This code could cause the problem of "cancellation of significant digits"...
	    COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs;
	  }else{
	    // Add the new index.
	    COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs;
	    for(int d = 0; d < DIM_STATE; d++){
	      STATE[num_of_terms][d] = working_state[d];
	    }
	    num_of_terms += 1;
	    if(num_of_terms >= MAX_ID){
	      fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
	      exit(1);
	    }
	  }
	}
      }
    }

    // Check the convergence.
    for(int i = 0; i < num_of_terms; i++){
      // Check whether the state is the target one or not.
      int is_target_state = TRUE;
      for(int d = 0; d < DIM_STATE; d++){
	if(STATE[i][d] != TARGET_STATE[d]){
	  is_target_state = FALSE;
	  break;
	}
      }
      if(is_target_state == TRUE){
	seq[0][M] = first_rate * COEFF[i] * final_rate;
	double coeff, intercept;
	// Evaluate the intercept (iteratively).
	for(int i_extra = 0; i_extra < num_diff-1; i_extra++){
	  coeff = (seq[i_extra][M] - seq[i_extra][M-1])/(1.0/(double)M - 1.0/(double)(M-1));
	  intercept = seq[i_extra][M] - coeff*(1.0/(double)M);
	  seq[i_extra+1][M] = intercept;
	}
	// Check whether a sequence of differentials is converged or not.
	// NOTE: The original sequence is excluded.
	for(int i_extra = 1; i_extra < num_diff; i_extra++){
	  char str1[256], str2[256];
	  sprintf(str1, "%+.10e", seq[i_extra][M]);
	  sprintf(str2, "%+.10e", seq[i_extra][M-1]);
	  int is_same = TRUE;
	  for(int ch = 0; ch <= CONV_DICIMAL_PLACE+2; ch++){
	    if(str1[ch] != str2[ch]){
	      is_same = FALSE;
	    }
	  }
	  if(is_same == TRUE){
	    is_conv = TRUE;
	  }
	}
	//fprintf(fp, "%3d", M);
	fprintf(fp, "%3d %e", M, 1/(double)M);
	for(int i = 0; i < num_diff; i++){
	  fprintf(fp, " %+.10e", seq[i][M]);
	}
	fprintf(fp, "\n");
	break;
      }
    }
    M += 1;
  }
  M = M-1;

  // Evaluate the weighted sum.
  // NOTE: The original sequence is excluded.
  double *w;
  w = (double*)clear_alloc(sizeof(double), num_diff+1);
  double eps = 1.0e-5;
  int start = 1;
  for(int i = start; i < num_diff; i++){
    w[i] = 1.0/(fabs(seq[i][M] - seq[i][M-1]) + eps);
  }
  double sum = 0.0;
  for(int i = start; i < num_diff; i++){
    sum += w[i];
  }
  for(int i = start; i < num_diff; i++){
    w[i] = w[i]/sum;
  }
  double total_result = 0.0;
  for(int i = start; i < num_diff; i++){
    total_result += w[i]*seq[i][M];
  }
  
  fclose(fp);
  free(w);
  free(working_state);
  free_2d(seq);
  return total_result;
}



double occur_event(int event, int *state, double dt_M){
  // Make the next state caused by the event.
  int is_valid_next_state = TRUE;
  for(int d = 0; d < DIM_STATE; d++){
    // If the next state is in the negative region, it is not valid.    
    if(state[d] + EVENT[event][d] < 0){
      is_valid_next_state = FALSE;
      break;
    }
  }
  if(is_valid_next_state == FALSE){
    return 0.0;
  }
  // Calculate the denominator.
  double rate = 0.0;
  for(int j = 0; j < LIST_LEN_E2ORI[ID_DIAGONAL_EVENT]; j++){
    int pre_e = LIST_START_INDEX_E2ORI[ID_DIAGONAL_EVENT] + j;
    double state_dependent_rate = 1.0;
    for(int d = 0; d < DIM_STATE; d++){
      for(int fact = 0; fact < ORI_INDEX4RATE[pre_e][d]; fact++){
	state_dependent_rate = state_dependent_rate * (double)(state[d] - fact);
      }
    }
    rate = rate + state_dependent_rate * ORI_EVENT_CONST_RATE[pre_e];
  }
  double denom_inv;
  denom_inv = 1.0/(1.0 - dt_M*rate);
  if(event == ID_DIAGONAL_EVENT){
    return denom_inv;
  }else{
    denom_inv = denom_inv*denom_inv;
    rate = 0.0;
    for(int j = 0; j < LIST_LEN_E2ORI[event]; j++){
      int pre_e = LIST_START_INDEX_E2ORI[event] + j;
      double state_dependent_rate = 1.0;
      for(int d = 0; d < DIM_STATE; d++){
	for(int fact = 0; fact < ORI_INDEX4RATE[pre_e][d]; fact++){
	  state_dependent_rate = state_dependent_rate * (double)(state[d] - fact);
	}
      }
      rate = rate + state_dependent_rate * ORI_EVENT_CONST_RATE[pre_e];
    }
    double numer;
    numer = - dt_M*rate;
    for(int d = 0; d < DIM_STATE; d++){
      state[d] += EVENT[event][d];
    }
    return - numer * denom_inv;
  }  
  return rate;
}


// -----------------------------------------------------------
// Subroutines for file open and memory allocation
// -----------------------------------------------------------
FILE *f_open(const char *path, const char *mode){
  static FILE *pf_ret;
  pf_ret = fopen(path, mode);
  if( pf_ret == NULL ){
    fprintf( stderr, "can't open a file, %s.\n", path );
    exit( EXIT_FAILURE );
  }
  return pf_ret;
}

void *clear_alloc(size_t size, size_t nobj){
  void *pret;
  pret = calloc(nobj, size);
  if (pret == NULL) {
    fprintf(stderr,"memory allocation failure!");
    exit(0);
  }
  return pret;
}

void **clear_alloc_2d(size_t size, size_t row, size_t column){
  void **ppret = NULL;
  size_t i;
  ppret = (void **) clear_alloc(sizeof(void *), row);
  ppret[0] = (void *) clear_alloc(size, row * column);
  for (i = 0; i < row; i++) {
    ppret[i] = (void *) ((size_t) ppret[0] + i * size * column);
  }
  return ppret;
}

