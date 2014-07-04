
#include <mistral_sat.hpp>

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
 

double *sorting_array;
int compar(const void *a, const void *b)
{
  if (sorting_array[*(int*)a]==sorting_array[*(int*)b])
    return 0;
  else
    if (sorting_array[*(int*)a] < sorting_array[*(int*)b])
      return 1;
    else
      return -1;
}
void initSort(double *sa)
{
  sorting_array = sa;
}


using namespace Mistral;
using namespace std;


void Mistral::print_literal(ostream& o, Literal l, unsigned int start_from, bool dir)
{

  o << (SIGN(l) ? "b" : "~b") << UNSIGNED(l)+start_from ;
  // if(dir)
  //   o << " == " << (l%2) ; 
  // else
  //   o << " =/= " << (1-(l%2)) ; 
  // //o << (l%2 ? "+" : "-") << UNSIGNED(l)+1;
}

void Mistral::print_clause(ostream& o, Clause *cl, unsigned int start_from)
{
  Clause& clause = *cl;
  o //<< " " << cl 
    << "(";
  for(unsigned int i=0; i<clause.size-1; ++i) {
    print_literal(o,clause[i], start_from);
    o << " v ";
  }
  print_literal(o,clause[clause.size-1],start_from);
  o << ")";
}


void SatSolver::set_parameters(SolverParameters& p) {
  params = p;
  set_policy(params.restart_policy);
}


void SatSolver::set_policy( const int policy )
{
  if( policy == LUBY )
    restart_policy = new Luby(params.restart_base);
  else
    restart_policy = new Geometric(params.restart_base, params.restart_factor);
}

SatSolver::SatSolver() 
{
  restart_policy = NULL;
}

SatSolver::SatSolver(const char* filename) 
{
  restart_policy = NULL;
  parse_dimacs(filename);
}

SatSolver::~SatSolver() 
{
  for(unsigned int i=0; i<learnt.size; ++i)
    free(learnt[i]);
  for(unsigned int i=0; i<base.size; ++i)
    free(base[i]);
  delete restart_policy;
}

int SatSolver::solve()
{
  usrand(params.seed);
  stats.start_time = get_run_time();
  if(!restart_policy)
    set_policy( params.restart_policy );
  // if(params.verbosity>1)
  //   cout << endl 
  // 	 << "c  ==================================[ Mistral (Sat module) ]===================================" << endl
  // 	 << "c  |      SEARCH  STATISTICS           |                  PROBLEM STATISTICS                   |" << endl
  // 	 << "c  |  Conflicts |       Nodes  CPUTime |  Atoms  Clauses (Size)   Learnt (Size)   Added (Size) |" << endl
  // 	 << "c  =============================================================================================" << endl;
  int result = (unsigned int)(stats.base_avg_size)+1;
  while(1) {
    //print_decisions(std::cout);
    if(params.shuffle) shuffle();
    if(params.dynamic_value) params.value_selection = randint(6);

    result = iterative_search();
    if(restart_policy) restart_policy->reset(params.restart_limit);
    if(params.verbosity>1)
      cout << "c  | " << setw(10) << stats.num_failures
	   << " | " << setw(10) << stats.num_nodes
	   << " " << setw(9) << (get_run_time() - stats.start_time)
	   << " | " << setw(6) << (state.size - stats.literals - 1)
	   << " " << setw(8) << base.size
	   << " " << setw(5) << setprecision(3) << stats.base_avg_size     
	   << "  " << setw(8) << learnt.size
	   << " " << setw(5) << setprecision(3) << stats.learnt_avg_size ;
    forget();
    if(params.verbosity>1)
      cout << "  " << setw(7) << learnt.size
	   << " " << setw(5) << setprecision(3) << stats.learnt_avg_size 
	   << "  |" << endl;
    if(result == LIMITOUT)
      status = UNKNOWN;
    else break;
  } 
  
  cout << "c  =============================================================================================" << endl;
  if( result == SAT ) 
    cout << "c  |                                    SATISFIABLE!                                           |" << endl;
  else if( result == UNSAT ) 
    cout << "c  |                                     UNSATISFIABLE!                                        |" << endl;
  else if( result == UNKNOWN ) 
    cout << "c  |                                       UNKNOWN!                                            |" << endl;
  double cputime = (get_run_time() - stats.start_time);
  if( cputime < 0.001 ) cputime = 0.001;
  cout << "c  =============================================================================================" << endl ;
  if(params.verbosity>0)    
    cout << "c  |      Conflicts  |          Nodes  |    CPU Time  | Conflicts/s |     Nodes/s | node/confl |" << endl
	 << "c  | " << setw(14) << stats.num_failures 
	 << "  | " << setw(14) << stats.num_nodes 
	 << "  | " << setw(11) << cputime
	 << "  | " << setw(11) << int(double(stats.num_failures)/cputime) 
	 << " | " << setw(11) << int(double(stats.num_nodes)/cputime) 
	 << " | " << setw(10) << (double(stats.num_nodes)/double(stats.num_failures))
	 << " |" << endl
	 << "c  =============================================================================================" << endl << endl;
  if(params.verbosity > 3) {
    if(result == SAT)
      {
	cout << "s SATISFIABLE\nv";
	for(unsigned int i=1; i<state.size; ++i) 
	  cout << " " << SIGN(state[i]) ;
	cout << endl;
      } else cout << "s UNSATISFIABLE" << endl;
    cout << "d ASSIGNMENTS " << stats.num_nodes << endl
	 << "d FAILS " << stats.num_failures << endl
      //<< "d BACKTRACKS " << stats.num_failures << endl
	 << "d NODES/s " << int(double(stats.num_nodes)/cputime) << endl;
  }


//   if( result == SAT )
//     print_decisions(cout);

  return status;
}

std::ostream& SatSolver::display(std::ostream& o) const {
  print_clauses(o);
  return o;
}

void SatSolver::print_all(ostream& o) const
{
  print_decisions( o );
  print_clauses( o );
  print_watchers( o );
}

void SatSolver::print_watchers(ostream& o, int beg, int end) const
{
  if( beg == NOVAL ) beg = -state.size;
  if( end == NOVAL ) end = state.size;
  
  Literal l;
  for(int i=beg; i<=end; ++i)
    {
      if(i) {
	l = (i<0 ? -2*(i+1) : (2*(i-1)+1));
	print_literal(o, l);
	o << " is watched by";
	o.flush();
	for(unsigned int j=0; j<is_watched_by[l].size; ++j) {
	  assert(is_watched_by[l][j]);
	  print_clause(cout, is_watched_by[l][j]);
	}
	o << endl;
      }
    }
}

void SatSolver::print_decisions(ostream& o, bool mode) const
{
  if( mode ) {
    cout << "assumptions: ";
    for(unsigned int i=0; i<assumptions.size; ++i) {
      o << " " ;
      Literal l = ((2*assumptions[i]) | SIGN(state[assumptions[i]]));
      print_literal(o, l);
      
      if(LEVEL(state[assumptions[i]]) != i)
	o << " warning - inconsistent level: " << LEVEL(state[assumptions[i]]);
    }
    o << endl;
  } else {
    unsigned int tlvl = 0;
    unsigned int clvl = decisions[0];
    for(unsigned int i=0; i<assumptions.size; ++i) {
      if(tlvl<decisions.size && i >= clvl) clvl = decisions[++tlvl];
	
      o << LEVEL(state[assumptions[i]]) // lvl[assumptions[i]]
	<< "\t" << tlvl
	<< "\t";
      print_literal(o, (2*assumptions[i]) | SIGN(state[assumptions[i]]));
      //o << "\t" << SIGN(state[assumptions[i]]) // polarity[assumptions[i]]
      o << "\t";
      o.flush();
      if( tlvl && reason[assumptions[i]] )
	print_clause(o, reason[assumptions[i]]);
      else if(tlvl)
	cout << "decision";
      else
	cout << "data";
      cout << endl;
    }
  }
}

void SatSolver::print_clauses(ostream& o) const
{
  o << "base (" << base.size << "):" << endl;
  for(unsigned int i=0; i<base.size; ++i)
    {
      o << "c" << i ;
      print_clause( o, base[i] );
      o << endl;
    }
  o << endl;
  o << "learnt (" << learnt.size << "):" << endl;
  for(unsigned int i=0; i<learnt.size; ++i)
    {
      o << "c" << i ;
      print_clause( o, learnt[i] );
      o << endl;
    }
  o << endl;
}

void SatSolver::init_vars(const int n, const int m)
{
  usrand(12345);

  status = UNKNOWN;
  decisions.initialise(0,n);
  assumptions.initialise(0,n);

  state.initialise(0,n);
  reason.initialise(0,n);
  //decisions_level.init(0,n+1);
  //assumptions_size.init(0,n+1);

  is_watched_by.initialise(0,2*n);
  activity.initialise(0,2*n);

  base.initialise(0,m);
  learnt.initialise(0,m);
  
//   init_num_atoms = n;
//   init_num_base = m;

  next_deduction = 0;

  visited.initialise(0, n-1, BitSet::empt);

  for(int i=0; i<n; ++i) {
    state.add(2*i + randint(2));
    reason[i] = NULL;
    
    assumptions[i] = i;
    decisions[i] = 0;
  }

  for(int i=0; i<2*n; ++i) {
    activity[i] = 0;
  }
}

void SatSolver::init_watchers()
{
  unsigned int i, j=base.size;

  for(i=0; i<2*state.size; ++i) {
    is_watched_by[i].initialise(0,(int)(2*activity[i]));
  }

  for(i=0; i<j; ++i) {
    Clause& clause = *(base[i]);
    is_watched_by[clause[0]].add(base[i]);
    is_watched_by[clause[1]].add(base[i]);
  }
} 

void SatSolver::parse_dimacs(const char* filename) 
{
  unsigned int LARGENUMBER = 131072;
  ifstream infile( filename );
  char c=' ';
  string word;
  int N, M, l=0, lit, cs;

  // skip comments
  infile >> c;
  while( c != 'p' ) {
    infile.ignore( LARGENUMBER, '\n' );
    infile >> c;
  }

  infile >> word;
  assert( word == "cnf" );
  
  // get number of atoms and clauses
  infile >> N;
  infile >> M;

  //init(N, M);
  init_vars(N, M);
  learnt_clause.initialise(0, N);

  for(int i=0; i<M; ++i)
    {
      learnt_clause.clear();
      do {
	infile >> l;
	if(l!=0) {
	  if(l>0) lit = (l-1)*2+1;
	  else lit = (l+1)*-2;
	  learnt_clause.add(lit);
	}
      } while(l && infile.good());


      if(params.init_activity == 1) {
	cs = learnt_clause.size;
	for(int j=cs-1; j>=0; --j) {
	  //std::cout << learnt_clause << " " << (params.activity_increment/cs) << std::endl;
	  activity[learnt_clause[j]] += (params.activity_increment/cs);
	}
      }

      add_clause( learnt_clause );
      if(params.checked) add_original_clause( learnt_clause );
    }

  init_watchers();

  if(params.normalize_activity != 0)
    normalize_activity(params.normalize_activity);
}

int SatSolver::check_solution()
{
  if(params.checked) {
    Atom x;
    bool correct = true;
    
    for(unsigned int i=0; correct && i<original.size; ++i) {
      bool satisfied = false;
      Clause& clause = *(original[i]);
      
      for(unsigned int j=0; !satisfied && j<clause.size; ++j) {
	x = ATOM(clause[j]);
	satisfied = (SIGN(x) == SIGN(clause[j]));
      }
      if(!satisfied) {
	print_clause(cerr, original[i]);
	cerr << endl;
      }
      
      correct = satisfied;
    }
    if(!correct) {
      cerr << "/!\\ The solution is not correct /!\\" << endl; 
      exit(1);
      return UNKNOWN;
    }
  }
  return SAT;
}


/*
Mistral::ConstraintClauseBase::ConstraintClauseBase(Vector< Variable >& scp) 
  : GlobalConstraint(scp) { 
  conflict = NULL;
  priority = 1;
}
*/
Mistral::ConstraintClauseBase::ConstraintClauseBase(Vector< Variable >& scp, bool __fd_variables, int st_from)
  : GlobalConstraint(scp), fd_variables(__fd_variables), start_from(st_from) {
  conflict = NULL;
  priority = 1;

  init_var_size = scp.size;


}

void Mistral::ConstraintClauseBase::mark_domain() {
  for(unsigned int i=0; i<scope.size; ++i) {
    ((Solver*)solver)->mark_non_convex(scope[i].id());
  }
}

void Mistral::ConstraintClauseBase::initialise() {
  //solver->base = this;

	ConstraintImplementation::initialise();

  for(unsigned int i=0; i<scope.size; ++i) {
    if(scope[i].is_bool())
      trigger_on(_VALUE_, scope[i]);
    else {
     // trigger_on(_NEVER_, scope[i]);
       std::cout << scope[i] << " " << scope[i].get_domain() << std::endl;
       exit(1);
    }
  }
  //set_idempotent(true);

  GlobalConstraint::initialise();

  is_watched_by.initialise(0,2*scope.size);
#ifdef _IMPROVE_UP
  clauses_of_literal.initialise(0,2*scope.size);
#endif
  // lit_activity.initialise(0,2*scope.size);
  // var_activity.initialise(0,scope.size);


  //reason_for = new Explanation*[scope.size];

  reason_for.initialise(0,scope.size);
  
  //Only when not lazy generation
  if (fd_variables && get_solver()->parameters.lazy_generation ){
	  std::cout << " c will enforce_nfc1 to false" << std::endl;
	  enforce_nfc1 = false;
  }

  will_be_kept.initialise(0,2000000,BitSet::empt);

  // for(unsigned int i=0; i<scope.size; ++i) {
  //   reason.add(NULL);
  // //   is_watched_by[i] = new Vector< Clause* >;
  // //   is_watched_by[i]-=scope[i].get_min();
  // }

  //reason = solver->reason.stack_;

}

Mistral::ConstraintClauseBase::~ConstraintClauseBase() {
  for(unsigned int i=0; i<clauses.size; ++i) {
    free(clauses[i]);
  }
  for(unsigned int i=0; i<learnt.size; ++i) {
    free(learnt[i]);
  }
}

void Mistral::ConstraintClauseBase::add(Variable x) {
  unsigned int idx = x.id();
  if(idx == scope.size) {
    scope.add(x);
    reason_for.add(NULL);

  } else if(idx > scope.size) {
    while(scope.capacity <= idx)
      scope.extendStack();
    scope[idx] = x;
    
    while(reason_for.capacity <= idx)
      reason_for.extendStack();
    reason_for[idx] = NULL;

    // while(is_watched_by.capacity <= 2*idx)
    //   is_watched_by.extendStack();
    // while(lit_activity.capacity <= 2*idx)
    //   lit_activity.extendStack();
    // while(var_activity.capacity <= idx)
    //   var_activity.extendStack();
  }

  while(is_watched_by.capacity <= 2*idx)
    is_watched_by.extendStack();

#ifdef _IMPROVE_UP
  while(clauses_of_literal.capacity <= 2*idx)
	  clauses_of_literal.extendStack();
#endif
  // while(lit_activity.capacity <= 2*idx)
  //   lit_activity.extendStack();
  // while(var_activity.capacity <= idx)
  //   var_activity.extendStack();
  
  //reason = solver->reason.stack_;
}




void Mistral::ConstraintClauseBase::extend_scope(Variable x){

	//	scope.add(x);

	//add(x);
	/*
		if (is_watched_by.capacity < (2*scope.size))
			is_watched_by.extendStack();
		if (reason_for.capacity < (scope.size))
			reason_for.extendStack();
	 */
	unsigned int idx = x.id()-start_from;
	if(idx == scope.size) {
		scope.add(x);
		//   _scope.add(x);
		//   reason_for.add(NULL);

		while(reason_for.capacity <= idx)
			reason_for.extendStack();
		reason_for[idx] = NULL;

	} else {
		//if(idx > scope.size) {

		std::cout << "idx > scope.size) ? " << std::endl;
		exit(1);

		while(scope.capacity <= idx)
			scope.extendStack();
		scope[idx] = x;

		while(_scope.capacity <= idx)
			_scope.extendStack();
		_scope[idx] = x;

		while(reason_for.capacity <= idx)
			reason_for.extendStack();
		reason_for[idx] = NULL;

		// while(is_watched_by.capacity <= 2*idx)
		//   is_watched_by.extendStack();
		// while(lit_activity.capacity <= 2*idx)
		//   lit_activity.extendStack();
		// while(var_activity.capacity <= idx)
		//   var_activity.extendStack();
	}

	while(is_watched_by.capacity <= 2*idx)
		is_watched_by.extendStack();

#ifdef _IMPROVE_UP
	while(clauses_of_literal.capacity <= 2*idx)
		clauses_of_literal.extendStack();
#endif
	// while(lit_activity.capacity <= 2*idx)
	//   lit_activity.extendStack();
	// while(var_activity.capacity <= idx)
	//   var_activity.extendStack();

	//reason = solver->reason.stack_;


	++_currentsize;
	if (_currentsize>= _capacity)
	{
		extend_vectors(5000);
		std::cout << " c _currentsize>= _capacity in ConstraintClauseBase" << std::endl;

	/*
	Event * tmp_event_type = new Event[scope.size];
	int * tmp_solution = new int[scope.size];
	Constraint* tmpself = new Constraint[scope.size];
	int*  tmpindex = new int[scope.size];

	for(unsigned int i=0; i<on.size; ++i)
		tmp_event_type[i] = event_type[i];

	delete [] event_type;
	event_type = tmp_event_type;


	for(unsigned int i=0; i<on.size; ++i)
		tmp_solution[i] = solution[i];

	delete [] solution;
	solution = tmp_solution;



	for(unsigned int i=0; i<on.size; ++i)
		tmpself[i] = self[i];

	delete [] self;
	self = tmpself;

	for(unsigned int i=0; i<on.size; ++i)
		tmpindex[i] = index[i];

	delete [] index;
	index = tmpindex;
	*/
	}

	event_type[scope.size-1] = NO_EVENT;
	solution [scope.size-1]  = 0;

	/*
	std:: cout << " \n BEFORE : scope.size " << scope.size << std::endl;
	std:: cout << "changes.size " << changes.size << std::endl;
	std:: cout << "changes.list_capacity " << changes.list_capacity << std::endl;
	std:: cout << "changes.index_capacity " << changes.index_capacity << std::endl;
	std:: cout << "changes " << changes << std::endl;
	 */
	if ((changes.list_capacity < scope.size) || (changes.index_capacity < scope.size)){
		changes.extend_lists();
	/*
	std:: cout << " \n After : scope.size " << scope.size << std::endl;
	std:: cout << "changes.size " << changes.size << std::endl;
	std:: cout << "changes.list_capacity " << changes.list_capacity << std::endl;
	std:: cout << "changes.index_capacity " << changes.index_capacity << std::endl;
	std:: cout << "changes " << changes << std::endl;
	 */



	//events.size = changes.size;
	events.index_capacity = changes.index_capacity;
	events.list_capacity = changes.list_capacity;
	events.list_ = changes.list_;
	events.index_ = changes.index_;
	}

	//	events.start_ = changes.start_;

	if (events.start_){
		std::cout << "start_ " << events.start_ << std::endl;
		std::cout << "exit start_ " << *events.start_ << std::endl;

		exit(1);
	}

	int i = scope.size -1;
	self[i] = Constraint(this, i|type);
	trigger_on(_VALUE_, scope[scope.size -1]);
	index[i] = on[i]->post(self[i]);
	consolidate_var(scope.size -1);
}





void Mistral::ConstraintClauseBase::add( Vector < Literal >& clause, double activity_increment ) {
 if(clause.size > 1) {
   Clause *cl = (Clause*)(Clause::Array_new(clause));
   clauses.add( cl );
   is_watched_by[clause[0]].add(cl);
   is_watched_by[clause[1]].add(cl);

#ifdef _IMPROVE_UP
   unsigned _size = clause.size;
   for (unsigned int i = 0; i< _size; ++i){
	   clauses_of_literal[clause[i]].add(cl);
   }
#endif

   // // should we split the increment?
   // activity_increment /= clause.size;

   // if(activity_increment > 0.0) {
   //   int i=clause.size;
   //   while(i--) {

   //     //std::cout << clause << " " << activity_increment << std::endl;

   //     lit_activity[clause[i]] += activity_increment;
   //     var_activity[UNSIGNED(clause[i])] += activity_increment;
   //   }
   // }

 } else {
   scope[UNSIGNED(clause[0])].set_domain(SIGN(clause[0]));
 }
}

void Mistral::ConstraintClauseBase::learn( Vector < Literal >& clause, double activity_increment, bool forget_immediately) {
 if(clause.size > 1) {
   Clause *cl = (Clause*)(Clause::Array_new(clause));
   learnt.add( cl );

   //std::cout << " c LEARN " << forget_immediately  << std::endl;
   // // should we split the increment?
   // activity_increment /= clause.size;

   // if(activity_increment > 0.0) {
   //   int i=clause.size;
   //   while(i--) {
   //     lit_activity[clause[i]] += activity_increment;
   //     var_activity[UNSIGNED(clause[i])] += activity_increment;
   //   }
   // }
   if (forget_immediately){
	   will_be_forgotten.add(learnt.size -1);
   }
   else{
	   is_watched_by[clause[0]].add(cl);
	   is_watched_by[clause[1]].add(cl);

#ifdef _IMPROVE_UP
	   unsigned _size = clause.size;
	   for (unsigned int i = 0; i< _size; ++i){
		   clauses_of_literal[clause[i]].add(cl);
	   }
#endif

   }

 } else {
   scope[UNSIGNED(clause[0])].set_domain(SIGN(clause[0]));
 }
}

void Mistral::ConstraintClauseBase::initialise_activity(double *lact, double *vact, double norm) {
  
  int i=clauses.size, j;

  // there are 2^n assignments in a clause, only one that falsifies it.
  double activity_increment;


  while(i--) {
    Clause& clause = *(clauses[i]);
    
    //activity_increment = norm / (1 << (clause.size-1));
    activity_increment = norm / (clause.size-1);
    activity_increment = norm / (clause.size-1);
    activity_increment = norm / (clause.size-1);
    j=clause.size;
    while(j--) {
      lact[NOT(clause[j])] += activity_increment;
      vact[UNSIGNED(clause[j])] += activity_increment;
    }
  }

}

int Mistral::ConstraintClauseBase::check( const int* sol ) const {
  unsigned int i, j;
  bool falsified=false;

  for(i=0; !falsified && i<clauses.size; ++i) {
    falsified = true;
    Clause& clause = *(clauses[i]);

    // std::cout << clause << " (" ;
    // for(j=0; j<clause.size; ++j) {
    //   //std::cout << " " << scope[UNSIGNED(clause[j])]
    //   std::cout << " ";
    //   print_literal(std::cout, clause[j]);
    // }
    // std::cout << " )" << std::endl;
    // for(j=0; j<clause.size; ++j) {
    //   std::cout << " " << (sol[UNSIGNED(clause[j])]) ;
    // }
    // std::cout << std::endl;

    for(j=0; falsified && j<clause.size; ++j) {
      falsified = //clause[j].check(sol[j]);
	(sol[UNSIGNED(clause[j])] != (int)SIGN(clause[j]));
    }

    //if(falsified) exit(1);
  }
  
  return falsified;
}

//#define _DEBUG_UNITPROP true

void Mistral::ConstraintClauseBase::start_over() {

	//std::cout << "init_var_size" << init_var_size <<  std::endl;

	on.size = init_var_size;
	reason_for.size =init_var_size ;

	//is_watched_by.
	for (int i =  init_var_size*2 ; i < ((2* scope.size) -1) ; ++i)
	  if (is_watched_by[i].size)
	  {
			std::cout << "ERROR  (is_watched_by[i].size) " <<  (is_watched_by[i].size) <<  std::endl;
			exit(1);
	  }


#ifdef _IMPROVE_UP
	for (int i =  init_var_size*2 ; i < ((2* scope.size) -1) ; ++i)
		if (clauses_of_literal[i].size)
		{
			std::cout << "ERROR  (clauses_of_literal[i].size) " <<  (clauses_of_literal[i].size) <<  std::endl;
			exit(1);
		}
#endif

	scope.size = init_var_size;
	_scope.size = init_var_size;

	for(int i=0; i<changes.index_capacity; ++i)
	{
		changes.index_[i] = i;
		if(i < (int)changes.list_capacity) changes.list_[i] = i;
	}

	changes.size=0;
	events.size = changes.size;
	events.index_capacity = changes.index_capacity;
	events.list_capacity = changes.list_capacity;
	events.list_ = changes.list_;
	events.index_ = changes.index_;
	_currentsize = init_var_size;

}

//#define _CHECKED_CLAUSES true

Mistral::PropagationOutcome Mistral::ConstraintClauseBase::propagate() {
  conflict=NULL;
  PropagationOutcome wiped = CONSISTENT;

//  return wiped;
  int x, v, cw;
  Literal p;

  while( !conflict && !changes.empty() ) {
    x = changes.pop();
    v = scope[x].get_min();
    p = NOT(2*x+v);

#ifdef _DEBUG_UNITPROP
    for(unsigned int i=0; i<solver->level; ++i) std::cout << " " ;
    std::cout << "propagate " ;
    print_literal(std::cout, NOT(p), start_from);
  //  std::cout << " that is : index :  " << x << " value : " <<v << " literal " << p;
    std::cout << " " << is_watched_by[p].size << std::endl;
    std::cout << " " << is_watched_by[NOT(p)].size << std::endl;
#endif


#ifdef _IMPROVE_UP

	unsigned int _size = clauses_of_literal[NOT(p)].size;

	for (unsigned int i = 0; i< _size; ++i){
	    update_clauses_of_literal(NOT(p), i);
	}

#endif

    cw = is_watched_by[p].size;
    while(cw-- && !conflict) {
      conflict = update_watcher(cw, p, wiped);
    }

/*
    cw = is_watched_by[NOT(p)].size;
    while(cw-- && !conflict) {
      conflict = update_watcher(cw, NOT(p), wiped);
    }
*/

  }

#ifdef _CHECKED_CLAUSES
  unsigned int i, j;
  if(IS_OK(wiped)) {
    bool violated = false;
    int num_literals = 0;
    for(i=0; !violated && num_literals != 1 && i<clauses.size; ++i) {
      Clause& clause = *(clauses[i]);
      violated = true;
      num_literals = clause.size;
      for(j=0; j<clause.size; ++j) {
	Atom a = UNSIGNED(clause[j]);
	if(scope[a].is_ground()) {
	  --num_literals;
	  violated &= !scope[a].equal(SIGN(clause[j]));
	} else violated = false;
      }
    }
    if(violated || num_literals==1) {
    	--i;
      std::cout << "unit propagation was not complete!! \n before printing, should we decrese i? i.e. --i? " << std::endl;
      print_clause(std::cout, clauses[i],start_from);
      std::cout << std::endl;
      exit(1);
    } else {

    	for(i=0; !violated && num_literals <= 1 && i<learnt.size; ++i) {
    		//A test if will_be_forgotten doesn't contain i
    		if (!will_be_forgotten.size  ||
    		(!will_be_forgotten.index_of_elt(i) &&  will_be_forgotten[0] != i)){

    		Clause& clause = *(learnt[i]);
    		violated = false;
    		num_literals = clause.size;
    		// if(clause.size == 23) {
    		//   std::cout << "check ";
    		//   print_clause(std::cout, learnt[i]);
    		//   std::cout << std::endl;
    		//   for(j=0; j<clause.size; ++j) {
    		//     Atom a = UNSIGNED(clause[j]);
    		//     std::cout << scope[a].get_domain() << " ";
    		//   }
    		//   std::cout << std::endl;
    		// }

    		for(j=0; j<clause.size; ++j) {
    			Atom a = UNSIGNED(clause[j]);
    			if(scope[a].is_ground()) {
    				--num_literals;
    				//	violated &= !scope[a].equal(SIGN(clause[j]));
    				if (scope[a].equal(SIGN(clause[j]))){
    					violated = false;
    					break;
    				}
    				else
    					violated= true;
    			} else violated = false;
    		}
    	}
    	}

    	if(violated  && num_literals<=1) {
    		--i;
    		std::cout << "unit propagation was not complete!!" << std::endl;
    		std::cout << "violated !" << violated <<std::endl;
    		std::cout << "num_literals " << num_literals <<std::endl;
    		//		std::cout << "level : " << get_solver()->level << std::endl;

    		print_clause(std::cout, learnt[i], start_from);
    		Clause& clause = *(learnt[i]);
    		std::cout << "\n Literals " <<  std::endl;

    		for (j = 0; j <clause.size ; ++j )
    			std::cout << "  " << clause[j] ;

    		std::cout << " \n  variables " <<  std::endl;

    		for (j = 0; j <clause.size; ++j )
    			//std::cout << get_solver()->get_id_boolean_variable(learnt[i][j])
    			std::cout << " " << scope[UNSIGNED(clause[j])].get_domain();


    		std::cout << std::endl;

    		std::cout << std::endl;



    		std::cout << " \n  printing watshers : " <<  std::endl;
    		std::cout << " \n learnt.size : " << learnt.size <<  std::endl;

    		//	for(i=0; i< is_watched_by.size ; ++i)
    		for(i=0; i< scope.size*2 ; ++i)
    		{
    			//	  Literal l =

    			if (is_watched_by[i].size)
    			{
    				std::cout << " literal : " << i << std::endl;;
    				print_literal(std::cout, i,start_from);
    				std::cout << " \n \n \n \n is watching :" << std::endl;;
    			}
    			//std::cout.flush();
    			for(unsigned int j=0; j<is_watched_by[i].size; ++j) {
    				// assert(is_watched_by[l][j]);
    				print_clause(cout, is_watched_by[i][j], start_from);
    				std::cout << "" << std::endl;;
    			}

    		}
    		std::cout << " END" << std::endl;;
    		exit(1);
    	}




    	/*
      for(i=0; !violated && num_literals != 1 && i<learnt.size; ++i) {
	Clause& clause = *(learnt[i]);
	violated = true;
	num_literals = clause.size;
	// if(clause.size == 23) {
	//   std::cout << "check ";
	//   print_clause(std::cout, learnt[i]);
	//   std::cout << std::endl;
	//   for(j=0; j<clause.size; ++j) {
	//     Atom a = UNSIGNED(clause[j]);
	//     std::cout << scope[a].get_domain() << " ";
	//   }
	//   std::cout << std::endl;
	// }

	for(j=0; j<clause.size; ++j) {
	  Atom a = UNSIGNED(clause[j]);
	  if(scope[a].is_ground()) {
	    --num_literals;
	    violated &= !scope[a].equal(SIGN(clause[j]));
	  } else violated = false;
	}
      }
      if(violated || ((violated) && num_literals==1)) {
    	  --i;
  		std::cout << "unit propagation was not complete!!" << std::endl;
		std::cout << "violated !" << violated <<std::endl;
		std::cout << "num_literals " << num_literals <<std::endl;
    	    //		std::cout << "level : " << get_solver()->level << std::endl;

	print_clause(std::cout, learnt[i], start_from);
	Clause& clause = *(learnt[i]);
	std::cout << "\n Literals " <<  std::endl;

	for (j = 0; j <clause.size ; ++j )
		std::cout << "  " << clause[j] ;

	std::cout << " \n  variables " <<  std::endl;

	for (j = 0; j <clause.size; ++j )
		//std::cout << get_solver()->get_id_boolean_variable(learnt[i][j])
		std::cout << " " << scope[UNSIGNED(clause[j])].get_domain();


	std::cout << std::endl;

	std::cout << std::endl;



	std::cout << " \n  printing watshers : " <<  std::endl;
	std::cout << " \n learnt.size : " << learnt.size <<  std::endl;

//	for(i=0; i< is_watched_by.size ; ++i)
	for(i=0; i< scope.size*2 ; ++i)
	{
		//	  Literal l =

		if (is_watched_by[i].size)
		{
			std::cout << " literal : " << i << std::endl;;
			print_literal(std::cout, i,start_from);
			std::cout << " \n \n \n \n is watching :" << std::endl;;
		}
		//std::cout.flush();
		for(unsigned int j=0; j<is_watched_by[i].size; ++j) {
			// assert(is_watched_by[l][j]);
			print_clause(cout, is_watched_by[i][j], start_from);
			std::cout << "" << std::endl;;
		}

	}
	std::cout << " END" << std::endl;;
	exit(1);
      }

    	 */
    }
  }
#endif


#ifdef _VERIFY_BEHAVIOUR_WHEN_LEARNING
  if (wiped != CONSISTENT)
  {
	  //  	  std::cout << "wiped == " << wiped << std::endl;
	  ( (Solver*) solver) ->__failure = this;
#ifdef _DEBUG_FAIL
	  std::cout << " fail : " << *this << std::endl;
#endif
  }
#endif

  if (conflict){
	  get_solver()-> update_failure_scope(conflict);
  }

  return wiped;
}


Mistral::Clause* Mistral::ConstraintClauseBase::update_watcher(const int cw, 
							       const Literal p,
							       PropagationOutcome& po)
{
  Clause *cl = is_watched_by[p][cw];
  Clause& clause = *cl;
  unsigned int j;

  Literal q, r;
  Variable v, w;
  int vb, wb;

#ifdef _DEBUG_WATCH
  std::cout << "update watchers for " << clause 
	    << " because " << (SIGN(p) ? "" : "~") << UNSIGNED(p)
	    << " <-> " << scope[UNSIGNED(p)] << " in " 
	    << scope[UNSIGNED(p)].get_domain() << std::endl;
#endif

  //ensure that p is the second watched lit
  if( clause[1] != p ) {
    q = clause[1];
    clause[0] = q;
    clause[1] = p;
  } else q = clause[0];
  v=scope[UNSIGNED(q)];
  vb=*(v.bool_domain);

  //check if the other watched lit is assigned
  //if( !v.is_ground() || v.get_min() != (int)SIGN(q) ) {
  if( vb==3 || vb>>1 != (int)SIGN(q) ) {

#ifdef _DEBUG_WATCH    
    std::cout << "  the second watcher does not satisfy the clause, we need a replacement" << std::endl;
#endif

    for(j=2; j<clause.size; ++j) {
      // for each literal r of the clause,
      r = clause[j];
      w = scope[UNSIGNED(r)];
      wb = *(w.bool_domain);

#ifdef _DEBUG_WATCH
      std::cout << "    what about " << (SIGN(r) ? "" : "~") << UNSIGNED(r)
		<< " <-> " << w << " in " << w.get_domain() << std::endl; 
#endif

      if( wb == 3 ) { // this literal is not set
	//if( !w.is_ground() ) { // this literal is not set
	// then it is a good candidate to replace p

	clause[1] = r;
	clause[j] = p;
	is_watched_by[p].remove(cw);
	is_watched_by[r].add(cl);

#ifdef _DEBUG_WATCH
	std::cout << "    ok!" // << clause << " " << (cl)
		  << std::endl;
#endif

	break;	
      }
      // if it is set true, then the clause is satisfied
      //else if( w.get_min() == (int)SIGN(r) ) {
      else if( wb>>1 == (int)SIGN(r) ) {

#ifdef _DEBUG_WATCH
	std::cout << "    ok! (satisfied)" << std::endl;
#endif

	clause[1] = r;
	clause[j] = p;
	is_watched_by[p].remove(cw);
	is_watched_by[r].add(cl);

	break;
      }
    }
      
    if( j == clause.size ) // no replacement could be found
      { 

#ifdef _DEBUG_WATCH
	std::cout << "  couldn't find a replacement!" << std::endl;
#endif

	//if( !v.is_ground() ) {
	if( vb == 3 ) {
	  // the last literal (other watched lit) is not set yet, we set it
	  //add_lit(q);

	  //std::cout << "    -> b" << UNSIGNED(q) << " = " << SIGN(q) << std::endl;

	  changes.add(UNSIGNED(q));
	  v.set_domain(SIGN(q));
	  //reason[UNSIGNED(q)] = cl;
	  //EXPL
	  reason_for[UNSIGNED(q)] = cl;
	  

#ifdef _DEBUG_UNITPROP
	  std::cout << "    -> " << v << " in " << v.get_domain() << std::endl;
#endif

	} else 
	  // it is set to false already, we fail
	  //if( v.get_min() != (int)SIGN(q) ) {
		//To do : we do not need this test!
	  if( vb>>1 != (int)SIGN(q) ) {

#ifdef _DEBUG_WATCH
	    std::cout << "    -> fail!" << std::endl;
#endif
	    po = FAILURE(UNSIGNED(q));

	    return cl;
	  }
      }
  }

  return NULL;
}

#ifdef _IMPROVE_UP
void Mistral::ConstraintClauseBase::update_clauses_of_literal(const Literal p, unsigned int i)
{
	//	unsigned int _size = clauses_of_literal[p].size;

	Clause *	cl = clauses_of_literal[p][i];
	Clause& clause = *cl;
	//Clause& clause = *cl;
	Literal q ;
	//Variable v;
	int vb;
	unsigned int position =0 ,  _s = clause.size;

	//for (unsigned int i = 0; i< _size; ++i){
	//cl = clauses_of_literal[p][i];
	//	Clause& clause = *cl;

	//		position =0 ;

	//		_s = clause.size;
	for (unsigned int j=0; j<_s; ++j)
		if (clause[j]==p){
			position=j;
			break;
		}

	//		if ((position==0) || (position ==1))
	//			return;
	if (position>1){

		q = clause[0];

		vb=*(scope[UNSIGNED(q)].bool_domain);

		//check if the other watched lit is assigned
		//if( !v.is_ground() || v.get_min() != (int)SIGN(q) ) {
		if( vb==3 || vb>>1 != (int)SIGN(q) ) {
			q = clause[1];
			//v=scope[UNSIGNED(q)];
			vb=*(scope[UNSIGNED(q)].bool_domain);
			if( vb==3 || vb>>1 != (int)SIGN(q) ) {
				//Here non of the first literals satisfy the clause!
				q = clause[0];
				//v=scope[UNSIGNED(q)];
				//vb=*(v.bool_domain);
				// this literal is not set
				//if( !w.is_ground() ) { // this literal is not set
				// then it is a good candidate to replace p

				clause[0] = p;
				clause[position] = q;
				is_watched_by[p].add(cl);
				is_watched_by[q].remove_elt(cl);
			}
		}
	}
	//}
}
#endif

void Mistral::ConstraintClauseBase::remove( const int cidx , bool static_forget)
{
  Clause *clause = learnt[cidx];

  // print_clause(std::cout, clause);
  // std::cout << std::endl;

  if (!static_forget){
	  is_watched_by[clause->data[0]].remove_elt( clause );
	  is_watched_by[clause->data[1]].remove_elt( clause );

#ifdef _IMPROVE_UP
	  unsigned _size = clause->size;
	  Clause& c = *clause;
	  for (unsigned int i = 0; i< _size; ++i){
		  clauses_of_literal[c[i]].remove_elt(clause);
	  }
#endif
  }

  //We suppose that we never remove a clause that will_be_kept!

  if (will_be_kept.fast_contain(learnt.size -1)){
	 // std::cout << "will_be_kept" << will_be_kept << std::endl;
	  //std::cout << "will_be_kept fast_contain " << learnt.size -1 << std::endl;


	  // here we update the index of learnt.size -1 since it will be kept!
	  will_be_kept.fast_remove(learnt.size -1);
	  will_be_kept.fast_add(cidx);
  }

  //std::cout << "  \n WILL REMOVE  \n "  << std::endl;
  //std::cout << " "  << learnt[cidx] << std::endl;

  learnt.remove( cidx );

  free(clause);
}


void Mistral::ConstraintClauseBase::static_forget(){

	/*std::cout << "\n \n \n \n  c static_forget " << std::endl;
	std::cout << " c will_be_forgotten size" << will_be_forgotten.size << std::endl;
	std::cout << " c learnt size" << learnt.size << std::endl;

	std::cout << " c will_be_forgotten " << will_be_forgotten << std::endl;
	std::cout << " c learnt " << learnt << std::endl;


	std::cout << " static forget " ;*/
	for (int i = (will_be_forgotten.size -1); i >=0 ; --i){
		remove(will_be_forgotten[i], true);
	}
	will_be_forgotten.clear();

	/*std::cout << " c will_be_forgotten " << will_be_forgotten << std::endl;
	std::cout << " c learnt " << learnt << std::endl;
	exit(1);
	 */

	//	std::cout << " c END hard_forget" << std::endl;
	//	std::cout << " c will_be_forgotten size" << will_be_forgotten.size << std::endl;
	//	std::cout << " c learnt size" << learnt.size << std::endl;

}


//#define _DEBUG_FORGET true


int Mistral::ConstraintClauseBase::forget(const double forgetfulness, 
					  const double *var_activity,
					  const double *lit_activity
					  //const Vector< double >& lit_activity
					  )
{

	//std::cout << "\n \n c start forget with " << learnt.size << std::endl;
	//TODO update remove with static forget
  int removed = 0;
  int * solution = get_solver()->last_solution_lb.stack_;


  int nlearnt = learnt.size;
  int keep=1;
  int i=0;
  Atom a;


  if( forgetfulness > 0.0 ) {
    // int nlearnt = learnt.size;
    double sa[nlearnt];
    Clause *tmp[nlearnt];
    int j, order[nlearnt], real_size;
    initSort(&(sa[0]));
    for(i=0; i<nlearnt; ++i)
      {
	order[i] = i;
	if(lit_activity) {
	  sa[i] = 0.0;
	  Clause& clause = *(learnt[i]);
	  real_size = j = clause.size;
	  while(j--) // THE ACTIVITY OF A LITERAL IS A MEASURE OF HOW MUCH IT IS "WANTED" BY THE FORMULA - SHORT CLAUSE WITH UNWANTED LITERALS ARE THEREFORE GOOD
	    {
	      a = UNSIGNED(clause[j]) + start_from;
	      // real_size is the number of free literal in hte clause.
	      // For correcteness, we need to not forget any clause that currently explains a literal.
	      // It seems like a good idea to keep clauses with small real size anyway (they matter the most right now)
	      //if(scope[UNSIGNED(clause[j])].is_ground()) --real_size;
	      if(solution[a] != -INFTY && solution[a] != (int)(SIGN(clause[j]))) --real_size;
	      //else 
	      sa[i] += (var_activity[a] + lit_activity[NOT(clause[j]) + (2*start_from)]);
	    }
	  // if(real_size) {
	  //   sa[i] /= (double)(real_size);
	  //   sa[i] /= (double)(clause.size *clause.size);
	  // //(real_size * real_size); //
	  // //(clause.size * clause.size);
	  // }
	  // else
	  //   sa[i] = INFTY;
	  sa[i] /= (double)((real_size+1) *clause.size *clause.size);
	  //sa[i] /= (double)(clause.size *clause.size *clause.size);
	} else {
	  sa[i] = 1.0/(double)(learnt[i]->size);
	}
      }

    keep = (int)((double)nlearnt * (1.0-forgetfulness));

    if(lit_activity)
      qsort(order, nlearnt, sizeof(int), compar);
    else {
      int swap;
      for(i=0; i<keep; ++i) {
	j = i + randint(nlearnt-i);
	swap = order[i];
	order[i] = order[j];
	order[j] = swap;
      }
    }

    for(i=0; i<nlearnt; ++i)
      tmp[i] = learnt[order[i]]; 
    for(i=0; i<nlearnt; ++i) {
      learnt[i] = tmp[i];
    }
    

#ifdef _DEBUG_FORGET
    //bool weird = true;
    std::cout << " BEFORE _DEBUG_FORGET, check if weight should include  var_activity \n" << std::endl;
    exit(1);


    for(i=nlearnt; --i>=0;) {
      double weight = 0;
      Clause& clause = *(learnt[i]);
      real_size = j = clause.size;
      while(j--)
	{
	  a = UNSIGNED(clause[j])+start_from;
	  if(solution[a] != -INFTY && solution[a] != (int)(SIGN(clause[j]))) --real_size;
	  // else 
	  weight += lit_activity[NOT(clause[j]) + (2*start_from)];
	}

      std::cout << setw(3) << i << ": " << learnt[i]->size << " " << real_size << " " << sa[order[i]] << "/" << weight/(double)((real_size+1) * clause.size * clause.size) << std::endl;
      if(i==keep) std::cout << "=================================\n";

      // if(weird && i<keep && i && sa[order[i]] != sa[order[i-1]]) {
      // 	std::cout << "=================================\n";
      // 	weird = false;
      // }
    }
#endif


    for(i=nlearnt; i>keep && sa[order[i-1]] != INFTY;) {
      removed += learnt[i-1]->size;
      remove( --i );
    }

    
    // /// PIECE OF CODE THAT I CAN'T UNDERSTAND!!!!
    // while(i>1 && sa[order[i-1]] != INFTY) {
    //   --i;
    //   if(sa[order[i]] == sa[order[i-1]]) {
	
    // 	// std::cout << "REMOVE  ";
    // 	// print_clause(std::cout, learnt[i]) ;
    // 	// std::cout << "\nBECAUSE ";
    // 	// print_clause(std::cout, learnt[i-1]) ;
    // 	// std::cout << std::endl;

    // 	avg_rem2_size += learnt[i]->size;

    // 	removed += learnt[i]->size;
    // 	remove( i );
    //   }
    // }

    // for(int kk=i; kk; ) {

    //   // std::cout << sa[order[kk-1]] << " keep " ;
    //   // print_clause(std::cout, learnt[kk-1]) ;
    //   // std::cout << std::endl;

    //   --kk;
      
    //   avg_kept_size += learnt[kk]->size;
    // }

  }

  //TODO improve static forget with parameters
  //Static forget
  nlearnt = learnt.size;
  if (nlearnt>12000){
	  std::cout << " c size limit in learntClauses : " << learnt.size << std::endl;
	  //removed= 0;
	 // nlearnt -=1000;
	  for(i=nlearnt; i>1 ;) {
		  removed += learnt[i-1]->size;
		  remove( --i );
	  }
  }

  //  std::cout << " c end forget with " << learnt.size << std::endl;
  //  std::cout << " c learnt" << learnt << std::endl;

  return removed;

}


std::ostream& Mistral::ConstraintClauseBase::display(std::ostream& os) const {
  os << " (";
  if(clauses.size>0) {
    if(clauses.size<100) {
      print_clause(os, clauses[0]);
      
      for(unsigned int i=1; i<clauses.size; ++i) {
	os << " " ;
	print_clause(os, clauses[i]);
      }
      

      // os << clauses[0];
      // for(unsigned int i=1; i<clauses.size; ++i)
      //   os << " " << clauses[i]  ;
    } else {
      std::cout << "many clauses";
    }
  }
  os << ")";
  //os << "nogoods";
  return os;
}






// Mistral::ConstraintExtClauseBase::ConstraintExtClauseBase(Vector< Variable >& scp) 
//   : GlobalConstraint(scp) { 
//   conflict = NULL;
//   priority = 1;
// }

// void Mistral::ConstraintExtClauseBase::mark_domain() {
//   for(unsigned int i=0; i<scope.size; ++i) {
//     ((Solver*)solver)->mark_non_convex(scope[i].id());
//   }
// }

// void Mistral::ConstraintExtClauseBase::initialise() {
//   //solver->base = this;

//   for(unsigned int i=0; i<scope.size; ++i) {
//     trigger_on(_VALUE_, scope[i]);
//   }
//   //set_idempotent(true);

//   GlobalConstraint::initialise();

//   is_watched_by.initialise(0,2*scope.size);
//   lit_activity.initialise(0,2*scope.size);
//   var_activity.initialise(0,scope.size);
//   //reason.initialise(0,scope.size);
  

//   // for(unsigned int i=0; i<scope.size; ++i) {
//   //   reason.add(NULL);
//   // //   is_watched_by[i] = new Vector< ExtClause* >;
//   // //   is_watched_by[i]-=scope[i].get_min();
//   // }

//   //reason = solver->reason.stack_;

// }

// Mistral::ConstraintExtClauseBase::~ConstraintExtClauseBase() {
//   for(unsigned int i=0; i<clauses.size; ++i) {
//     free(clauses[i]);
//   }
//   for(unsigned int i=0; i<learnt.size; ++i) {
//     free(learnt[i]);
//   }
// }

// void Mistral::ConstraintExtClauseBase::add(Variable x) {
//   unsigned int idx = x.id();
//   if(idx == scope.size) {
//     scope.add(x);
//     //reason.add(NULL);

//   } else if(idx > scope.size) {
//     while(scope.capacity <= idx)
//       scope.extendStack();
//     scope[idx] = x;

//     // // while(reason.capacity <= idx)
//     // //   reason.extendStack();
//     // // reason[idx] = NULL;

//     // while(is_watched_by.capacity <= 2*idx)
//     //   is_watched_by.extendStack();
//     // while(lit_activity.capacity <= 2*idx)
//     //   lit_activity.extendStack();
//     // while(var_activity.capacity <= idx)
//     //   var_activity.extendStack();
//   }

//     while(is_watched_by.capacity <= 2*idx)
//       is_watched_by.extendStack();
//     while(lit_activity.capacity <= 2*idx)
//       lit_activity.extendStack();
//     while(var_activity.capacity <= idx)
//       var_activity.extendStack();

//   //reason = solver->reason.stack_;
// }

// void Mistral::ConstraintExtClauseBase::add( Vector < Literal >& clause, double activity_increment ) {
//  if(clause.size > 1) {
//    ExtClause *cl = (ExtClause*)(ExtClause::Array_new(clause));
//    clauses.add( cl );
//    is_watched_by[clause[0]].add(cl);
//    is_watched_by[clause[1]].add(cl);

//    // should we split the increment?
//    activity_increment /= clause.size;

//    if(activity_increment > 0.0) {
//      int i=clause.size;
//      while(i--) {

//        //std::cout << clause << " " << activity_increment << std::endl;

//        lit_activity[clause[i]] += activity_increment;
//        var_activity[UNSIGNED(clause[i])] += activity_increment;
//      }
//    }

//  } else {
//    scope[UNSIGNED(clause[0])].set_domain(SIGN(clause[0]));
//  }
// }

// void Mistral::ConstraintExtClauseBase::learn( Vector < Literal >& clause, double activity_increment ) {
//  if(clause.size > 1) {
//    ExtClause *cl = (ExtClause*)(ExtClause::Array_new(clause));
//    learnt.add( cl );

//    // // should we split the increment?
//    // activity_increment /= clause.size;

//    // if(activity_increment > 0.0) {
//    //   int i=clause.size;
//    //   while(i--) {
//    //     lit_activity[clause[i]] += activity_increment;
//    //     var_activity[UNSIGNED(clause[i])] += activity_increment;
//    //   }
//    // }

//    is_watched_by[clause[0]].add(cl);
//    is_watched_by[clause[1]].add(cl);
//  } else {
//    scope[UNSIGNED(clause[0])].set_domain(SIGN(clause[0]));
//  }
// }

// int Mistral::ConstraintExtClauseBase::check( const int* sol ) const {
//   unsigned int i, j;
//   bool falsified=false;

//   for(i=0; !falsified && i<clauses.size; ++i) {
//     falsified = true;
//     ExtClause& clause = *(clauses[i]);

//     // std::cout << clause << " (" ;
//     // for(j=0; j<clause.size; ++j) {
//     //   //std::cout << " " << scope[UNSIGNED(clause[j])]
//     //   std::cout << " ";
//     //   print_literal(std::cout, clause[j]);
//     // }
//     // std::cout << " )" << std::endl;
//     // for(j=0; j<clause.size; ++j) {
//     //   std::cout << " " << (sol[UNSIGNED(clause[j])]) ;
//     // }
//     // std::cout << std::endl;

//     for(j=0; falsified && j<clause.size; ++j) {
//       falsified = //clause[j].check(sol[j]);
// 	(sol[UNSIGNED(clause[j])] != (int)SIGN(clause[j]));
//     }

//     if(falsified) exit(1);
//   }
  
//   return falsified;
// }

// //#define _DEBUG_UNITPROP true



// ///! THIS VERSION NO LONGER WORKS BECAUSE TH ARRAY SEQUENCE IS NOT UPDATED IMMEDIATLY ON CHANGE ////
// // Mistral::PropagationOutcome Mistral::ConstraintExtClauseBase::propagate() {
// // //   conflict=NULL;
// // //   PropagationOutcome wiped = CONSISTENT;

// // //   int x, v, cw;
// // //   Literal p;

// // //   while( !conflict && !changes.empty() ) {
// // //     x = changes.pop();
// // //     v = scope[x].get_min();
// // //     p = NOT(2*x+v);

// // // #ifdef _DEBUG_UNITPROP
// // //     for(unsigned int i=0; i<solver->level; ++i) std::cout << " " ;
// // //     std::cout << "propagate " ;
// // //     print_literal(std::cout, NOT(p));
// // //     std::cout << " " << is_watched_by[p].size << std::endl;
// // // #endif

// // //     cw = is_watched_by[p].size;
// // //     while(cw-- && !conflict) {
// // //       conflict = update_watcher(cw, p, wiped);
// // //     }
// // //   }

// // //   return wiped;

// //   conflict=NULL;
// //   Variable *assumptions = ((Solver*)solver)->sequence.list_;
// //   int i=0, n=changes.size, index = ((Solver*)solver)->sequence.size;
// //   PropagationOutcome wiped = CONSISTENT;

// //   int x, v, cw;
// //   Literal p;

// //   //std::cout << solver->sequence << std::endl;
// //   //std::cout << "+++++++++++++++++" << std::endl;

// //   while( !conflict && i<n ) {
// //     x = changes[i];
// //     v = scope[x].get_min();
// //     p = NOT(2*x+v);

// // #ifdef _DEBUG_UNITPROP
// //     for(int j=0; j<solver->level; ++j) std::cout << " " ;
// //     std::cout << "propagate1 " ;
// //     print_literal(std::cout, NOT(p));
// //     std::cout << " " << is_watched_by[p].size << std::endl;
// // #endif

// //     cw = is_watched_by[p].size;
// //     while(cw-- && !conflict) {
// //       conflict = update_watcher(cw, p, wiped);
// //     }

// //     ++i;
// //   }


// //   //std::cout << "::::::::::::::::+" << std::endl;

// //   //std::cout << solver->sequence << std::endl;

// //   //std::cout << changes

// //   //n = solver->sequence.size;
// //   while( !conflict && --index>=(int)(((Solver*)solver)->sequence.size) ) {

// //     //std::cout << index << std::endl;

// //     x = assumptions[index].id();
// //     v = assumptions[index].get_min();
// //     p = NOT(2*x+v);

// // #ifdef _DEBUG_UNITPROP
// //     for(int j=0; j<solver->level; ++j) std::cout << " " ;
// //     std::cout << "propagate2 " ;
// //     print_literal(std::cout, NOT(p));
// //     std::cout << " " << is_watched_by[p].size << std::endl;
// // #endif

// //     cw = is_watched_by[p].size;
// //     while(cw-- && !conflict) {
// //       conflict = update_watcher(cw, p, wiped);
// //     }
// //   }

// //   //std::cout << "]]]]]]]]]]]]]]]]+" << std::endl;

// //   return wiped;
// // }


// Mistral::PropagationOutcome Mistral::ConstraintExtClauseBase::propagate() {
//   conflict=NULL;
//   PropagationOutcome wiped = CONSISTENT;

//   int x, v, cw;
//   Literal p;

//   while( !conflict && !changes.empty() ) {
//     x = changes.pop();
//     v = scope[x].get_min();
//     p = NOT(2*x+v);

// #ifdef _DEBUG_UNITPROP
//     for(unsigned int i=0; i<solver->level; ++i) std::cout << " " ;
//     std::cout << "propagate " ;
//     print_literal(std::cout, NOT(p));
//     std::cout << " " << is_watched_by[p].size << std::endl;
// #endif

//     cw = is_watched_by[p].size;
//     while(cw-- && !conflict) {
//       conflict = update_watcher(cw, p, wiped);
//     }
//   }

//   return wiped;

// //   conflict=NULL;
// //   Variable *assumptions = ((Solver*)solver)->sequence.list_;
// //   int i=0, n=changes.size, index = ((Solver*)solver)->sequence.size;
// //   PropagationOutcome wiped = CONSISTENT;

// //   int x, v, cw;
// //   Literal p;

// //   //std::cout << solver->sequence << std::endl;
// //   //std::cout << "+++++++++++++++++" << std::endl;

// //   while( !conflict && i<n ) {
// //     x = changes[i];
// //     v = scope[x].get_min();
// //     p = NOT(2*x+v);

// // #ifdef _DEBUG_UNITPROP
// //     for(int j=0; j<solver->level; ++j) std::cout << " " ;
// //     std::cout << "propagate1 " ;
// //     print_literal(std::cout, NOT(p));
// //     std::cout << " " << is_watched_by[p].size << std::endl;
// // #endif

// //     cw = is_watched_by[p].size;
// //     while(cw-- && !conflict) {
// //       conflict = update_watcher(cw, p, wiped);
// //     }

// //     ++i;
// //   }


// //   //std::cout << "::::::::::::::::+" << std::endl;

// //   //std::cout << solver->sequence << std::endl;

// //   //std::cout << changes

// //   //n = solver->sequence.size;
// //   while( !conflict && --index>=(int)(((Solver*)solver)->sequence.size) ) {

// //     //std::cout << index << std::endl;

// //     x = assumptions[index].id();
// //     v = assumptions[index].get_min();
// //     p = NOT(2*x+v);

// // #ifdef _DEBUG_UNITPROP
// //     for(int j=0; j<solver->level; ++j) std::cout << " " ;
// //     std::cout << "propagate2 " ;
// //     print_literal(std::cout, NOT(p));
// //     std::cout << " " << is_watched_by[p].size << std::endl;
// // #endif

// //     cw = is_watched_by[p].size;
// //     while(cw-- && !conflict) {
// //       conflict = update_watcher(cw, p, wiped);
// //     }
// //   }

// //   //std::cout << "]]]]]]]]]]]]]]]]+" << std::endl;

// //   return wiped;
// }


// Mistral::ExtClause* Mistral::ConstraintExtClauseBase::update_watcher(const int cw, 
// 							       const Literal p,
// 							       PropagationOutcome& po)
// {
//   ExtClause *cl = is_watched_by[p][cw];
//   ExtClause& clause = *cl;
//   unsigned int j;

//   Literal q, r;
//   Variable v, w;
//   int vb, wb;

// #ifdef _DEBUG_WATCH
//   std::cout << "update watchers for " << clause 
// 	    << " because " << (SIGN(p) ? "" : "~") << UNSIGNED(p)
// 	    << " <-> " << scope[UNSIGNED(p)] << " in " 
// 	    << scope[UNSIGNED(p)].get_domain() << std::endl;
// #endif

//   //ensure that p is the second watched lit
//   if( clause[1] != p ) {
//     q = clause[1];
//     clause[0] = q;
//     clause[1] = p;
//   } else q = clause[0];
//   v=scope[UNSIGNED(q)];
//   vb=*(v.bool_domain);

//   //check if the other watched lit is assigned
//   //if( !v.is_ground() || v.get_min() != (int)SIGN(q) ) {
//   if( vb==3 || vb>>1 != (int)SIGN(q) ) {

// #ifdef _DEBUG_WATCH    
//     std::cout << "  the second watcher does not satisfy the clause, we need a replacement" << std::endl;
// #endif

//     for(j=2; j<clause.size; ++j) {
//       // for each literal r of the clause,
//       r = clause[j];
//       w = scope[UNSIGNED(r)];
//       wb = *(w.bool_domain);

// #ifdef _DEBUG_WATCH
//       std::cout << "    what about " << (SIGN(r) ? "" : "~") << UNSIGNED(r)
// 		<< " <-> " << w << " in " << w.get_domain() << std::endl; 
// #endif

//       if( wb == 3 ) { // this literal is not set
// 	//if( !w.is_ground() ) { // this literal is not set
// 	// then it is a good candidate to replace p

// 	clause[1] = r;
// 	clause[j] = p;
// 	is_watched_by[p].remove(cw);
// 	is_watched_by[r].add(cl);

// #ifdef _DEBUG_WATCH
// 	std::cout << "    ok!" // << clause << " " << (cl)
// 		  << std::endl;
// #endif

// 	break;	
//       }
//       // if it is set true, then the clause is satisfied
//       //else if( w.get_min() == (int)SIGN(r) ) {
//       else if( wb>>1 == (int)SIGN(r) ) {

// #ifdef _DEBUG_WATCH
// 	std::cout << "    ok! (satisfied)" << std::endl;
// #endif

// 	break;
//       }
//     }
      
//     if( j == clause.size ) // no replacement could be found
//       { 

// #ifdef _DEBUG_WATCH
// 	std::cout << "  couldn't find a replacement!" << std::endl;
// #endif

// 	//if( !v.is_ground() ) {
// 	if( vb == 3 ) {
// 	  // the last literal (other watched lit) is not set yet, we set it
// 	  //add_lit(q);

// 	  //std::cout << "    -> b" << UNSIGNED(q) << " = " << SIGN(q) << std::endl;

// 	  changes.add(UNSIGNED(q));
// 	  v.set_domain(SIGN(q));
// 	  //reason[UNSIGNED(q)] = cl;
// 	  //EXPL
// 	  reason_for[UNSIGNED(q)] = cl;
	  

// #ifdef _DEBUG_UNITPROP
// 	  std::cout << "    -> " << v << " in " << v.get_domain() << std::endl;
// #endif

// 	} else 
// 	  // it is set to false already, we fail
// 	  //if( v.get_min() != (int)SIGN(q) ) {
// 	  if( vb>>1 != (int)SIGN(q) ) {

// #ifdef _DEBUG_WATCH
// 	    std::cout << "    -> fail!" << std::endl;
// #endif
// 	    po = FAILURE(UNSIGNED(q));

// 	    return cl;
// 	  }
//       }
//   }

//   return NULL;
// }

// void Mistral::ConstraintExtClauseBase::remove( const int cidx )
// {
//   ExtClause *clause = learnt[cidx];

//   is_watched_by[clause->data[0]].remove_elt( clause );
//   is_watched_by[clause->data[1]].remove_elt( clause );
//   learnt.remove( cidx );

//   free(clause);
// }


// void Mistral::ConstraintExtClauseBase::forget(double forgetfulness)
// {
//   if( forgetfulness > 0.0 ) {
//     int nlearnt = learnt.size;
//     double sa[nlearnt];
//     ExtClause *tmp[nlearnt];
//     int i, j, order[nlearnt];
//     initSort(&(sa[0]));
//     for(i=0; i<nlearnt; ++i)
//       {
// 	order[i] = i;
// 	sa[i] = 0.0;
// 	ExtClause& clause = *(learnt[i]);
// 	j=clause.size;
// 	while(j--) // THE ACTIVITY OF A LITERAL IS A MEASURE OF HOW MUCH IT IS "WANTED" BY THE FORMULA - SHORT CLAUSE WITH UNWANTED LITERALS ARE THEREFORE GOOD
// 	  //sa[i] += var_activity[UNSIGNED(clause[j])];
// 	  //sa[i] += lit_activity[clause[j]];
// 	  sa[i] += lit_activity[NOT(clause[j])];
// 	sa[i] /= (clause.size * clause.size);

// 	//std::cout << sa[i] << " ";
//       }
//     //std::cout << std::endl;
//     qsort(order, nlearnt, sizeof(int), compar);
//     for(i=0; i<nlearnt; ++i)
//       tmp[i] = learnt[order[i]]; 
//     for(i=0; i<nlearnt; ++i) {
//       learnt[i] = tmp[i];
//     }
    
//     int keep = (int)((double)nlearnt * (1.0-forgetfulness));
    
//     for(i=nlearnt; i>keep;)
//       remove( --i );
//     while(i>1) {
//       --i;
//       if(sa[order[i]] == sa[order[i-1]])
// 	remove( i );
//     }
//   }
// }


// std::ostream& Mistral::ConstraintExtClauseBase::display(std::ostream& os) const {
//   os << " (";
//   if(clauses.size>0) {
//     os << clauses[0];
//     for(unsigned int i=1; i<clauses.size; ++i)
//       os << " " << clauses[i]  ;
//   }
//   os << ")";
//   //os << "nogoods";
//   return os;
// }







// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ExtClause& x) {
//   os << "(" << (SIGN(x[0]) ? "" : "~") << UNSIGNED(x[0]) ;
//   for(unsigned int i=1; i<x.size; ++i) {
//     os << " " << (SIGN(x[i]) ? "" : "~") << UNSIGNED(x[i]) ;
//   }
//   os << ")";
//   return os;
// }

// std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ExtClause* x) {
//   return os << (*x);
// }


std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Clause& x) {
  os << "(" << (SIGN(x[0]) ? "" : "~") << UNSIGNED(x[0]) ;
  for(unsigned int i=1; i<x.size; ++i) {
    os << " " << (SIGN(x[i]) ? "" : "~") << UNSIGNED(x[i]) ;
  }
  os << ")";
  return os;
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Clause* x) {
  return os << (*x);
}

std::ostream& Mistral::operator<< (std::ostream& os, Mistral::SatSolver& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, Mistral::SatSolver* x) {
  return x->display(os);
}

