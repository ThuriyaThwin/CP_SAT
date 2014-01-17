
#include <mistral_search.hpp>
#include <mistral_variable.hpp>
#include <mistral_constraint.hpp>


// 111111111111111 <--> 32767 <--> 7FFF
#define  VARIABLE_MASK 0x7FFF;
// 11111111111111111000000000000000 <--> 4294934528 <--> FFFF8000
#define  VALUE_MASK 0xFFFF8000;

//#define ENCODE_BOUND_LITERAL(id_variable, value) (((value << 15) | id_variable);)
//#define GET_VARIABLE_FROM_LITERAL(literal) (literal &VARIABLE_MASK;)
//#define GET_VALUE_FROM_LITERAL(literal) (literal >> 15 ;)


using namespace std;
using namespace Mistral;

//with 32 bits :
//The first bit is the sign : 0 lower bound; 1 : upper bound.
//The next 16 bits for values : --> biggest value is 65535.
//the last 14 bits are for the variable_id : biggest variable_id is 32767.
/*
inline unsigned int encode_bound_literal (int id_variable, int value,int sign) {return ( (sign << 31) | (value << 15)   | id_variable);}
inline int get_variable_from_literal (unsigned int literal) { return ( 0x7FFF & literal) ;}
inline int get_value_from_literal (unsigned int literal) {return ( (literal & 0x7FFFFFFF) >> 15);}
inline int get_sign_from_literal (unsigned int literal) {return (literal >> 31);}
inline bool is_upper_bound (unsigned int literal) {return (literal >> 31) ;}
inline bool is_lower_bound (unsigned int literal) {return 1- (literal >> 31) ;}

inline bool is_a_bound_literal (int literal) {return (literal < 0) || (literal > 0x7FFF ) ;}
inline int get_id_boolean_variable (int literal, int start_from) {return (literal + start_from) ;}
 */

void testrewrite1()
{

	std::cout << "\n \n test rewrite equalities" << std::endl;

	Variable X(0,10);
	Variable Y(-5,5);
	Variable Z(-100,10);
	Variable A(0, 10);

	Solver s;
	//	s.parameters.verbosity =88;
	s.add(X==Y);

	std:: cout << s << std::endl;
	s.rewrite();
	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();
	//	cout << s.statistics << endl;
	if(result)
	{
		Solution solx( s.variables);
		//	cout << " X  solution: " << solx << endl;
	}

}

void testrewrite2()
{

	std::cout << "\n \n test rewrite with expressions" << std::endl;

	Variable X(0,10);
	Variable Y(-5,5);
	Variable Z(-100,10);
	Variable A(0, 10);

	Solver s;
	//	s.parameters.verbosity =88;
	s.add((X+Z)==(Y+A));

	std:: cout << s << std::endl;
	s.rewrite();
	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();
	//	cout << s.statistics << endl;
	if(result)
	{
		Solution solx( s.variables);
		//	cout << " X  solution: " << solx << endl;
	}

}




void testrewriteNotEq()
{

	std::cout << "\n \n test rewrite NotEq" << std::endl;

	Variable X(0,10);
	Variable Y(-5,5);
	Solver s;
	//	s.parameters.verbosity =88;
	s.add(X==Y);
	s.add(X!=Y);

	std:: cout << s << std::endl;
	s.rewrite();
	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();
	cout << s.statistics << endl;

}



void testrewriteEqANDNotEq()
{

	std::cout << "\n \n test rewrite equalities" << std::endl;

	Variable X(0,10);
	Variable Y(-5,5);
	Variable Z(-100,10);

	Solver s;
	//	s.parameters.verbosity =88;
	s.add(X==Y);
	s.add(Z==Y);
	s.add(X!=Z);
	//s.add(F!=A);

	std:: cout << s << std::endl;
	s.rewrite();
	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();
	cout << s.statistics << endl;

}



void testVariableRangeDomainHistory()
{

	std::cout << "\n \n ttestVariableRangeDomainHistory" << std::endl;
	Variable X(0,10);
	Variable Y(-5,5);
	Variable Z(-100,10);

	Solver s;
	//	s.parameters.verbosity =88;
	s.add(X+7 < Y);
	//	s.add(Z==Y);
	s.add(X!=Z);
	//s.add(F!=A);

	std:: cout << s << std::endl;
	//s.rewrite();
	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();
	cout << s.statistics << endl;
}


void testDisjunctive()
{

	std::cout << "\n \n testDisjunctive" << std::endl;
	Variable X(0,10221);
	Variable Y(0,5443);
	Variable Z(0,10992);
	VarArray disjunts;
	Solver s;
	//	s.parameters.verbosity =88;
	s.add(X+7 <= Y);
	//	s.add(Z==Y);
	//disjunts.add(Precedence(X,Z,8,82));
	//s.add(F!=A);

	//	Free(disjuncts[0]);

	s.add(disjunts[0]);

	std:: cout << s << std::endl;
	//s.rewrite();
	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();
	cout << s.statistics << endl;
}


int encodeLiteral(int variable_id, int value)
{
	//	cout << "var " << variable_id << endl;
	//	cout << "value " << value << endl;
	int literal = value << 15 ;

	//	std::cout << "value shifted " << literal << std::endl;
	literal |= variable_id;
	std::cout << "final shifted " << literal << std::endl;


	//cout << "value mask " <<(int) VALUE_MASK << std::endl;

	return (value << 15) | variable_id;
	//	return literal;
}

/*
int get_variable_from_literal(int literal)
{
	return literal &VARIABLE_MASK;
}

int get_value_from_literal(int literal)
{
	return literal >> 15;
}

 */

/*
void test_encode_literal()
{

	int var = 31730;
	int val = 65535;
	int sign = 1;

	std::cout << "var " << var<< std::endl;
	std::cout << "value " << val << std::endl;
	std::cout << "sign " << sign << std::endl;

	unsigned int lit =encode_bound_literal(var,val, sign) ;


	std::cout << "AGAIN \n " <<  std::endl;
	std::cout << "var " << var<< std::endl;
	std::cout << "value " << val << std::endl;
	std::cout << "sign " << sign << std::endl;

	//	int a = (l & 0x80000000 ) >> 31;
	//	std::cout << "sign of value " << a << std::endl;
	//	std::cout << "sign of int " << sizeof (int) << std::endl;

	std::cout << "literal " << lit << std::endl;

	std::cout << " \n \n reverse sence : \n get var  " << get_variable_from_literal(lit) << std::endl;
	std::cout << " get value  " << get_value_from_literal(lit) << std::endl;
	std::cout << " get sign  " << get_sign_from_literal(lit) << std::endl;

	if ( is_upper_bound(lit))
		std::cout << " is upper bound   " << is_upper_bound(lit) << std::endl;
	else
		std::cout << " is lower bound    " << is_lower_bound(lit) << std::endl;

	int a = -1038;
	std::cout << " \n \n \n int a =  " << a << std::endl;
	std::cout << " unsigned int a =  " << (unsigned int ) a << std::endl;
	/*
	std::cout << " signed   " << a << std::endl;
	lit = a;
	std::cout << " unsigned X 1  :  " << lit << std::endl;
	lit = (unsigned int) a;
	std::cout << " unsigned X 2 :  " << lit << std::endl;

	if ( is_a_bound_literal(a))
	{
		std::cout << " is a Bound literal ! s.t. :"  << std::endl;
		std::cout << "var_id  " << get_variable_from_literal((unsigned int) a) << std::endl;
		std::cout << " get value  " << get_value_from_literal((unsigned int) a) << std::endl;
		std::cout << " get sign  " << get_sign_from_literal((unsigned int) a) << std::endl;
		if ( is_upper_bound((unsigned int) a))
			std::cout << " is upper bound   " << is_upper_bound((unsigned int) a) << std::endl;
		else
			std::cout << " is lower bound    " << is_lower_bound((unsigned int) a) << std::endl;
	}
	else
	{
		std::cout << " is not a bound literal "<< std::endl;
		//	std::cout << " get_variable " << get_id_boolean_variable(a,1872) << std::endl;

	}

}
 */

void test_pointer()
{

	int *  aa = new int [5];
	int size = 5;


	int a = 3;
	int * t = new int [3];

	t[0]=0;
	t[1]=1;
	t[2]=2;

	while (a--)
	{
		std::cout << " a : " << a << std::endl;
	}
	std::cout << " AFTER 1 : " << a << std::endl;

	a=3;
	while (a--)
	{
		std::cout << " a : " << a << std::endl;
		std::cout << " t[a] : " << t[a] << std::endl;
	}
	std::cout << " AFTER 1 : " << a << std::endl;





	while (size--)
	{
		aa[size] = size;
		std::cout << " a : " << aa[size] << std::endl;
	}

	std::cout << " XXX "  << std::endl;

	int * b = & aa[4];
	std::cout << " b : " << b << std::endl;
	std::cout << " *b : " << *b << std::endl;
	*b = 87;

	std::cout << " b : " << b << std::endl;
	std::cout << " *b : " << *b << std::endl;

	//	b++* =
	size = 18;
	while (size--)
	{
		std::cout << " a : " << aa[size] << std::endl;
	}

}



class A
{
public:
	virtual void f(){cout << "A::f()" << endl;}
};

class B : public A
{
public:
	void f(){cout << "B::f()" << endl;}
};


void test_dynamic_cast()
{


	A a;
	B b;
	a.f();        // A::f()
	b.f();        // B::f()

	A *pA = &a;
	B *pB = &b;
	pA->f();      // A::f()
	pB->f();      // B::f()

	pA = &b;
	// pB = &a;      // not allowed
	if (dynamic_cast<B*>(&a)) // allowed but it returns NULL
	{
		std::cout << "yes? \n" ;

	}
	else
		std::cout << "NO \n" ;

}



void test_while()
{

	int a = 7;

	do
	{

		std::cout << " a : " << a << std::endl;
	} while( --a );

	std::cout << " AFTER END a : " << a << std::endl;

}



void test_while2()
{

	int a = 7;

	Vector<int > v;

	v.add(87);
	v.add(8);
	v.add(7);
	v.add(5);
	v.add(9);

	int idx = v.size;

	while (v[--idx]!= 19)
	{

	}

	std::cout << " AFTER Wile a : v[idx]< " << v[idx]<< std::endl;

	std::cout << " AFTER Wile a :idx  " << idx<< std::endl;

}




void makeACounterExample()
{

	std::cout << "\n \n makeACounterExample" << std::endl;

	Variable X1L(0,100);
	Variable X1U(0,100);
	Variable X2L(0,100);
	Variable X2U(0,100);
	Variable X3L(0,100);
	Variable X3U(0,100);

	Variable Y1L(0,100);
	Variable Y1U(0,100);
	Variable Y2L(0,100);
	Variable Y2U(0,100);
	Variable Y3L(0,100);
	Variable Y3U(0,100);

	Variable Z1L(0,100);
	Variable Z1U(0,100);
	Variable Z2L(0,100);
	Variable Z2U(0,100);
	Variable Z3L(0,100);
	Variable Z3U(0,100);

	Variable Kx1(1,100);
	Variable Kx2(1,100);
	Variable Kx3(1,100);

	Variable Ky1(1,100);
	Variable Ky2(1,100);
	Variable Ky3(1,100);

	Variable Kz1(1,100);
	Variable Kz2(1,100);
	Variable Kz3(1,100);


	//New variables :
	Variable X2U_(0,100);
	Variable X1L_(0,100);
	Variable Y1U_(0,100);
	Variable Y2L_(0,100);


	Variable Y3U_(0,100);
	Variable X3U_(0,100);

	Variable Z1L_(0,100);
	Variable Z2L_(0,100);
	Variable Z3L_(0,100);




	Solver s;
	//	s.parameters.verbosity =88;
	//	s.add(X<= (Y+2));
	//	s.add(Z <=(Y+1));
	//s.add(X!=Z);
	//s.add(F!=A);

	//Bounds constraints
	s.add(X1L <=X1U);
	s.add(X2L <=X2U);
	s.add(X3L <=X3U);

	s.add(Y1L <=Y1U);
	s.add(Y2L <=Y2U);
	s.add(Y3L <=Y3U);

	s.add(Z1L <=Z1U);
	s.add(Z2L <=Z2U);
	s.add(Z3L <=Z3U);


	//Sequencing constraints
	s.add(X2L >= (X1L + Kx1));
	s.add(X2U >= (X1U + Kx1));

	s.add(X3L >= (X2L + Kx2));
	s.add(X3U >= (X2U + Kx2));

	s.add(Y1L >= (Y2L + Ky2));
	s.add(Y1U >= (Y2U + Ky2));

	s.add(Y3L >= (Y1L + Ky1));
	s.add(Y3U >= (Y1U + Ky1));

	//On Z
	s.add(Z2L >= (Z1L + Kz1));
	s.add(Z2U >= (Z1U + Kz1));

	s.add(Z1L >= (Z3L + Kz3));
	s.add(Z1U >= (Z3U + Kz3));



	//Disjunctions

	//Decided ones :

	s.add(X1L >= (Z1L + Kz1));
	s.add(X1U >= (Z1U + Kz1));

	s.add(Y2L >= (Z2L + Kz2));
	s.add(Y2U >= (Z2U + Kz2));

	s.add(X3L >= (Y3L + Ky3));
	s.add(X3U >= (Y3U + Ky3));

	//Not decided yet

	s.add(X1U >= (Y1L + Ky1));
	s.add(Y1U >= (X1L + Kx1));

	s.add(X2U >= (Y2L + Ky2));
	s.add(Y2U >= (X2L + Kx2));

	s.add(Z3U >= (X3L + Kx3));
	s.add(X3U >= (Z3L + Kz3));

	s.add(Z3U >= (Y3L + Ky3));
	s.add(Y3U >= (Z3L + Kz3));


	s.add(Z2U >= (X2L + Kx2));
	s.add(X2U >= (Z2L + Kz2));

	s.add(Z1U >= (Y1L + Ky1));
	s.add(Y1U >= (Z1L + Kz1));


	// new
	s.add (X2U_ <= X2U);
	s.add (X1L_ >=  X1L);
	s.add (Y1U_ <= Y1U );
	s.add ( Y2L_ >= Y2L);

	s.add (Y3U_ <= Y3U );
	s.add (X3U_ <= X3U );

	s.add (Z1L_>=Z1L);
	s.add ( Z2L_>= Z2L);
	s.add ( Z3L_ >= Z3L);


	//domain consistency

	s.add (X2U_ >= X2L);
	s.add (X1L_ <=  X1U);
	s.add (Y1U_ >= Y1L );
	s.add ( Y2L_ <= Y2U);

	s.add (Y3U_ >= Y3L );
	s.add (X3U_ >= X3L );

	s.add (Z1L_ <=  Z1U);
	s.add ( Z2L_ <=  Z2U);
	s.add ( Z3L_ <= Z3U);


	// Propagated Sequencing Constraints

	//	s.add(Y3L_ >= (Y1L_ + Ky1));
	s.add(Y3U_ >= (Y1U_ + Ky1));
	//	s.add(Y3L_ >= (Y1L_ + Ky1));
	s.add(X3U_ >= (X2U_ + Kx2));


	s.add(X1L_ >= (Z1L_ + Kz1));
	s.add(Y2L_ >= (Z2L_ + Kz2));

	s.add(Z2L_ >= (Z1L_ + Kz1));

	s.add(Z1L_ >= (Z3L_ + Kz3));



	//Fail on x2> y2
	s.add (X2U_ < (Y2L_ + Ky2)) ;
	//Fail on (x_1 < y_1
	s.add (Y1U_ < (X1L_ + Kx2)) ;


	//Propagate the new decision

	s.add(Z3U >= (X3U_ + Kx3));
	s.add(Z3L_ >= (X3L + Kx3));



	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();

	//	cout << " X  result: " << result << endl;
	//std:: cout << s << std::endl;

	//	std:: cout << X.get_solution_int_value() << std::endl;
	//	std:: cout << Y.get_solution_int_value() << std::endl;
	//	std:: cout << Z.get_solution_int_value() << std::endl;

	if(result)
	{
		Solution solx( s.variables);
		//		cout << "  solution Found " << solx << endl;
		cout << "  solution Found " << endl;
		cout << endl;


		cout << " Kx :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << Kx1.get_solution_int_value() << std::endl;
		std:: cout << Kx2.get_solution_int_value() << std::endl;
		std:: cout << Kx3.get_solution_int_value() << std::endl;

		cout << endl;
		std:: cout << Ky1.get_solution_int_value() << std::endl;
		std:: cout << Ky2.get_solution_int_value() << std::endl;
		std:: cout << Ky3.get_solution_int_value() << std::endl;

		cout << endl;
		std:: cout << Kz1.get_solution_int_value() << std::endl;
		std:: cout << Kz2.get_solution_int_value() << std::endl;
		std:: cout << Kz3.get_solution_int_value() << std::endl;

		cout << "\n  X :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << X1L.get_solution_int_value() << std::endl;
		std:: cout << X1U.get_solution_int_value() << std::endl;

		std:: cout << X2L.get_solution_int_value() << std::endl;
		std:: cout << X2U.get_solution_int_value() << std::endl;
		std:: cout << X3L.get_solution_int_value() << std::endl;
		std:: cout << X3U.get_solution_int_value() << std::endl;

		cout << "\n Y :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << Y1L.get_solution_int_value() << std::endl;
		std:: cout << Y1U.get_solution_int_value() << std::endl;
		std:: cout << Y2L.get_solution_int_value() << std::endl;
		std:: cout << Y2U.get_solution_int_value() << std::endl;
		std:: cout << Y3L.get_solution_int_value() << std::endl;
		std:: cout << Y3U.get_solution_int_value() << std::endl;

		cout << " \n Z :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << Z1L.get_solution_int_value() << std::endl;
		std:: cout << Z1U.get_solution_int_value() << std::endl;
		std:: cout << Z2L.get_solution_int_value() << std::endl;
		std:: cout << Z2U.get_solution_int_value() << std::endl;
		std:: cout << Z3L.get_solution_int_value() << std::endl;
		std:: cout << Z3U.get_solution_int_value() << std::endl;




		cout << " \n The new values  :  "<< endl;
		cout << endl;
		cout << endl;

		std::cout << "X2U_" << X2U_.get_solution_int_value() << std::endl;
		std::cout << "X1L_" <<  X1L_.get_solution_int_value() << std::endl;
		std::cout << "Y1U_" <<  Y1U_.get_solution_int_value() << std::endl;
		std::cout << "Y2L_" <<  Y2L_.get_solution_int_value() << std::endl;

		cout << " \n less important :  "<< endl;
		cout << endl;
		cout << endl;

		std::cout << "Y3U_" <<  Y3U_.get_solution_int_value() << std::endl;
		std::cout << "X3U_" <<  X3U_.get_solution_int_value() << std::endl;
		std::cout << "Z1L_" <<  Z1L_.get_solution_int_value() << std::endl;
		std::cout << "Z2L_" <<  Z2L_.get_solution_int_value() << std::endl;
		std::cout << "Z3L_" <<  Z3L_.get_solution_int_value() << std::endl;


	}
	else
		cout << " UNSAT! " << endl;

	//	cout << s.statistics << endl;


}

//inline bool is_a_latest_upper_bound (Literal literal) {return 3 & (literal >> 30) ;}
//inline bool is_a_latest_lower_bound (Literal literal) {return 2 & (literal >> 30) ;}
//inline bool is_a_latest_bound_literal (Literal literal) {return  (literal >> 31) ;}
//id_var should not exceed 30 bits ..
//
//inline Literal encode_latest_bound_literal (unsigned int id_variable,int sign) {return ( ((2+sign) << 30) | id_variable);}
//inline int get_variable_from_latest_literal (Literal literal) { return ( 0x3FFFFFFF & literal) ;}

void test_encode_latest_literal()
{

	int var = 31234730;

	int sign = 0;

	std::cout << "var " << var<< std::endl;
	std::cout << "sign " << sign << std::endl;

	Literal lit =encode_latest_bound_literal(var, sign) ;



	std::cout << "literal " << lit << std::endl;

	std::cout << " \n \n reverse sence : \n get var  " << get_variable_from_latest_literal(lit) << std::endl;
	//std::cout << " get value  " << get_value_from_literal(lit) << std::endl;
	//std::cout << " get sign  " << get_sign_from_literal(lit) << std::endl;

	//	if ( is_a_latest_upper_bound(lit))
	std::cout << " is an upper bound   " << is_a_latest_upper_bound (lit) << std::endl;
	//	else
	std::cout << " is a lower bound    " << is_a_latest_lower_bound(lit) << std::endl;



}


void  check_BitSet (){


	BitSet visitedbound;
	BitSet visitedbound22;

	visitedbound.initialise(0,3,BitSet::empt);
	visitedbound22.initialise(0,1023);
	visitedbound22.extend(6);

	std::cout << "visited22 " << std::endl;
	std::cout << "visited 22" << 	visitedbound22 <<std::endl;
	std::cout << "visited 22 size " << 	visitedbound22.size() <<std::endl;

	visitedbound22.fast_add(3);

	std::cout << "AFTERRRR  " << std::endl;
	std::cout << "visited 22" << 	visitedbound22 <<std::endl;
	std::cout << "visited 22 size " << 	visitedbound22.size() <<std::endl;


	visitedbound22.clear();



	std::cout << "AFTERRRR  " << std::endl;
	std::cout << "visited 22" << 	visitedbound22 <<std::endl;
	std::cout << "visited 22 size " << 	visitedbound22.size() <<std::endl;


	std::cout << "visited " << std::endl;
	std::cout << "visited " << 	visitedbound <<std::endl;
	std::cout << "visited size " << 	visitedbound.size() <<std::endl;



	//	i =
	//	if (visitedbound.fast)
	visitedbound.fast_add(2);

	std::cout << "visited 2 ? " << 	visitedbound.fast_contain(2) <<std::endl;

	visitedbound.fast_add(1);
	std::cout << "visited size " << 	visitedbound.size() <<std::endl;
	std::cout << "visited  " << 	visitedbound <<std::endl;

	std::cout << "visited 0 ? " << 	visitedbound.fast_contain(0) <<std::endl;
	std::cout << "visited 1 ? " << 	visitedbound.fast_contain(1) <<std::endl;
	std::cout << "visited 2 ? " << 	visitedbound.fast_contain(2) <<std::endl;
	std::cout << "visited 3 ? " << 	visitedbound.fast_contain(3) <<std::endl;

	visitedbound.clear();
	std::cout << "After Clear size  " << 	visitedbound.size() <<std::endl;

	std::cout << "visited 0 ? " << 	visitedbound.fast_contain(0) <<std::endl;
	std::cout << "visited 1 ? " << 	visitedbound.fast_contain(1) <<std::endl;
	std::cout << "visited 2 ? " << 	visitedbound.fast_contain(2) <<std::endl;
	std::cout << "visited 3 ? " << 	visitedbound.fast_contain(3) <<std::endl;


}



void test_encode_latest_literal0()
{

	unsigned int var = 0;

	int sign = 1;

	std::cout << "var " << var<< std::endl;
	std::cout << "sign " << sign << std::endl;

	Literal lit =encode_latest_bound_literal(var, sign) ;



	std::cout << "literal " << lit << std::endl;

	std::cout << " \n \n reverse sence : \n get var  " << get_variable_from_latest_literal(lit) << std::endl;
	//std::cout << " get value  " << get_value_from_literal(lit) << std::endl;
	//std::cout << " get sign  " << get_sign_from_literal(lit) << std::endl;

	//	if ( is_a_latest_upper_bound(lit))
	std::cout << " is an upper bound   " << is_a_latest_upper_bound (lit) << std::endl;
	//	else
	std::cout << " is a lower bound    " << is_a_latest_lower_bound(lit) << std::endl;



}


//1 <= 2
void make_precedence(Solver *s, Variable *lb1, Variable *lb2, Variable *ub1, Variable *ub2,int k )
{

	s->add((*lb2) >= ((*lb1) + k));
	s->add((*ub2) >= ((*ub1) + k));


}
void make_precedence(Solver *s, Variable *lb1, Variable *lb2, Variable *ub1, Variable *ub2,Variable * k )
{

	s->add((*lb2) >= ((*lb1) + (*k)));
	s->add((*ub2) >= ((*ub1) + (*k)));


}
void make_disjunction(Solver *s, Variable *lb1,  Variable *lb2,  Variable *ub1, Variable *ub2,  int k1, int k2)
{

	s->add((*ub1) >= ((*lb2) + k2));
	s->add((*ub2) >= ((*lb1) + k1));

}




void make_disjunction(Solver *s, Variable *lb1,  Variable *lb2,  Variable *ub1, Variable *ub2,  Variable * k1, Variable * k2)
{

	s->add((*ub1) >= ((*lb2) + (*k2)));
	s->add((*ub2) >= ((*lb1) + (*k1)));

}



void make_disjunction(Solver *s, Variable *lb1,  Variable *lb2,  Variable *ub1, Variable *ub2,  Variable * k1, int k2)
{

	s->add((*ub1) >= ((*lb2) + k2));
	s->add((*ub2) >= ((*lb1) + (*k1)));

}


void make_disjunction(Solver *s, Variable *lb1,  Variable *lb2,  Variable *ub1, Variable *ub2,  int k1, Variable * k2)
{

	s->add((*ub1) >= ((*lb2) + (*k2)));
	s->add((*ub2) >= ((*lb1) + k1));

}
void make_Domain_Consistency(Solver *s, Variable *l, Variable *u)
{
	s->add((*l) <= (*u));
}

void make_lowerbound_propagation(Solver *s, Variable *newlb, Variable *oldlb)
{
	s->add((*newlb) >= (*oldlb));
}
void make_upperbound_propagation(Solver *s, Variable *newub, Variable *oldub)
{

	s->add((*newub) <= (*oldub));
}



void makeACounterExampleLearningwithClausesandrevisitingAdisjunction()
{
	/*

	std::cout << "\n \n makeACounterExample" << std::endl;

	Variable X1L(0,100);
	Variable X1U(0,100);
	Variable X2L(0,100);
	Variable X2U(0,100);
	Variable X3L(0,100);
	Variable X3U(0,100);

	Variable Y1L(0,100);
	Variable Y1U(0,100);
	Variable Y2L(0,100);
	Variable Y2U(0,100);
	Variable Y3L(0,100);
	Variable Y3U(0,100);

	Variable Z1L(0,100);
	Variable Z1U(0,100);
	Variable Z2L(0,100);
	Variable Z2U(0,100);
	Variable Z3L(0,100);
	Variable Z3U(0,100);

	Variable Kx1(1,100);
	Variable Kx2(1,100);
	Variable Kx3(1,100);

	Variable Ky1(1,100);
	Variable Ky2(1,100);
	Variable Ky3(1,100);

	Variable Kz1(1,100);
	Variable Kz2(1,100);
	Variable Kz3(1,100);


	//New variables :
	Variable X2U_(0,100);
	Variable X1L_(0,100);
	Variable Y1U_(0,100);
	Variable Y2L_(0,100);


	Variable Y3U_(0,100);
	Variable X3U_(0,100);

	Variable Z1L_(0,100);
	Variable Z2L_(0,100);
	Variable Z3L_(0,100);




	Solver s;
	//	s.parameters.verbosity =88;
	//	s.add(X<= (Y+2));
	//	s.add(Z <=(Y+1));
	//s.add(X!=Z);
	//s.add(F!=A);

	//Bounds constraints
	s.add(X1L <=X1U);
	s.add(X2L <=X2U);
	s.add(X3L <=X3U);

	s.add(Y1L <=Y1U);
	s.add(Y2L <=Y2U);
	s.add(Y3L <=Y3U);

	s.add(Z1L <=Z1U);
	s.add(Z2L <=Z2U);
	s.add(Z3L <=Z3U);


	//Sequencing constraints
	s.add(X2L >= (X1L + Kx1));
	s.add(X2U >= (X1U + Kx1));

	s.add(X3L >= (X2L + Kx2));
	s.add(X3U >= (X2U + Kx2));

	s.add(Y1L >= (Y2L + Ky2));
	s.add(Y1U >= (Y2U + Ky2));

	s.add(Y3L >= (Y1L + Ky1));
	s.add(Y3U >= (Y1U + Ky1));

	//On Z
	s.add(Z2L >= (Z1L + Kz1));
	s.add(Z2U >= (Z1U + Kz1));

	s.add(Z1L >= (Z3L + Kz3));
	s.add(Z1U >= (Z3U + Kz3));



	//Disjunctions

	//Decided ones :

	s.add(X1L >= (Z1L + Kz1));
	s.add(X1U >= (Z1U + Kz1));

	s.add(Y2L >= (Z2L + Kz2));
	s.add(Y2U >= (Z2U + Kz2));

	s.add(X3L >= (Y3L + Ky3));
	s.add(X3U >= (Y3U + Ky3));

	//Not decided yet

	s.add(X1U >= (Y1L + Ky1));
	s.add(Y1U >= (X1L + Kx1));

	s.add(X2U >= (Y2L + Ky2));
	s.add(Y2U >= (X2L + Kx2));

	s.add(Z3U >= (X3L + Kx3));
	s.add(X3U >= (Z3L + Kz3));

	s.add(Z3U >= (Y3L + Ky3));
	s.add(Y3U >= (Z3L + Kz3));


	s.add(Z2U >= (X2L + Kx2));
	s.add(X2U >= (Z2L + Kz2));

	s.add(Z1U >= (Y1L + Ky1));
	s.add(Y1U >= (Z1L + Kz1));


	// new
	s.add (X2U_ <= X2U);
	s.add (X1L_ >=  X1L);
	s.add (Y1U_ <= Y1U );
	s.add ( Y2L_ >= Y2L);

	s.add (Y3U_ <= Y3U );
	s.add (X3U_ <= X3U );

	s.add (Z1L_>=Z1L);
	s.add ( Z2L_>= Z2L);
	s.add ( Z3L_ >= Z3L);


	//domain consistency

	s.add (X2U_ >= X2L);
	s.add (X1L_ <=  X1U);
	s.add (Y1U_ >= Y1L );
	s.add ( Y2L_ <= Y2U);

	s.add (Y3U_ >= Y3L );
	s.add (X3U_ >= X3L );

	s.add (Z1L_ <=  Z1U);
	s.add ( Z2L_ <=  Z2U);
	s.add ( Z3L_ <= Z3U);


	// Propagated Sequencing Constraints

	//	s.add(Y3L_ >= (Y1L_ + Ky1));
	s.add(Y3U_ >= (Y1U_ + Ky1));
	//	s.add(Y3L_ >= (Y1L_ + Ky1));
	s.add(X3U_ >= (X2U_ + Kx2));


	s.add(X1L_ >= (Z1L_ + Kz1));
	s.add(Y2L_ >= (Z2L_ + Kz2));

	s.add(Z2L_ >= (Z1L_ + Kz1));

	s.add(Z1L_ >= (Z3L_ + Kz3));



	//Fail on x2> y2
	s.add (X2U_ < (Y2L_ + Ky2)) ;
	//Fail on (x_1 < y_1
	s.add (Y1U_ < (X1L_ + Kx2)) ;


	//Propagate the new decision

	s.add(Z3U >= (X3U_ + Kx3));
	s.add(Z3L_ >= (X3L + Kx3));



	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();

	//	cout << " X  result: " << result << endl;
	//std:: cout << s << std::endl;

	//	std:: cout << X.get_solution_int_value() << std::endl;
	//	std:: cout << Y.get_solution_int_value() << std::endl;
	//	std:: cout << Z.get_solution_int_value() << std::endl;

	if(result)
	{
		Solution solx( s.variables);
		//		cout << "  solution Found " << solx << endl;
		cout << "  solution Found " << endl;
		cout << endl;


		cout << " Kx :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << Kx1.get_solution_int_value() << std::endl;
		std:: cout << Kx2.get_solution_int_value() << std::endl;
		std:: cout << Kx3.get_solution_int_value() << std::endl;

		cout << endl;
		std:: cout << Ky1.get_solution_int_value() << std::endl;
		std:: cout << Ky2.get_solution_int_value() << std::endl;
		std:: cout << Ky3.get_solution_int_value() << std::endl;

		cout << endl;
		std:: cout << Kz1.get_solution_int_value() << std::endl;
		std:: cout << Kz2.get_solution_int_value() << std::endl;
		std:: cout << Kz3.get_solution_int_value() << std::endl;

		cout << "\n  X :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << X1L.get_solution_int_value() << std::endl;
		std:: cout << X1U.get_solution_int_value() << std::endl;

		std:: cout << X2L.get_solution_int_value() << std::endl;
		std:: cout << X2U.get_solution_int_value() << std::endl;
		std:: cout << X3L.get_solution_int_value() << std::endl;
		std:: cout << X3U.get_solution_int_value() << std::endl;

		cout << "\n Y :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << Y1L.get_solution_int_value() << std::endl;
		std:: cout << Y1U.get_solution_int_value() << std::endl;
		std:: cout << Y2L.get_solution_int_value() << std::endl;
		std:: cout << Y2U.get_solution_int_value() << std::endl;
		std:: cout << Y3L.get_solution_int_value() << std::endl;
		std:: cout << Y3U.get_solution_int_value() << std::endl;

		cout << " \n Z :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << Z1L.get_solution_int_value() << std::endl;
		std:: cout << Z1U.get_solution_int_value() << std::endl;
		std:: cout << Z2L.get_solution_int_value() << std::endl;
		std:: cout << Z2U.get_solution_int_value() << std::endl;
		std:: cout << Z3L.get_solution_int_value() << std::endl;
		std:: cout << Z3U.get_solution_int_value() << std::endl;




		cout << " \n The new values  :  "<< endl;
		cout << endl;
		cout << endl;

		std::cout << "X2U_" << X2U_.get_solution_int_value() << std::endl;
		std::cout << "X1L_" <<  X1L_.get_solution_int_value() << std::endl;
		std::cout << "Y1U_" <<  Y1U_.get_solution_int_value() << std::endl;
		std::cout << "Y2L_" <<  Y2L_.get_solution_int_value() << std::endl;

		cout << " \n less important :  "<< endl;
		cout << endl;
		cout << endl;

		std::cout << "Y3U_" <<  Y3U_.get_solution_int_value() << std::endl;
		std::cout << "X3U_" <<  X3U_.get_solution_int_value() << std::endl;
		std::cout << "Z1L_" <<  Z1L_.get_solution_int_value() << std::endl;
		std::cout << "Z2L_" <<  Z2L_.get_solution_int_value() << std::endl;
		std::cout << "Z3L_" <<  Z3L_.get_solution_int_value() << std::endl;


	}
	else
		cout << " UNSAT! " << endl;

	//	cout << s.statistics << endl;


	 */

	Variable X1LB_(0,100);
	Variable X1UB_(0,100);

	Variable X2LB_(0,100);
	Variable X2UB_(0,100);

	Variable X3LB_(0,100);
	Variable X3UB_(0,100);

	Variable X4LB_(0,100);
	Variable X4UB_(0,100);

	Variable X5LB_(0,100);
	Variable X5UB_(0,100);


	Variable Y1LB_(0,100);
	Variable Y1UB_(0,100);

	Variable Y2LB_(0,100);
	Variable Y2UB_(0,100);

	Variable Y3LB_(0,100);
	Variable Y3UB_(0,100);

	Variable Y4LB_(0,100);
	Variable Y4UB_(0,100);

	Variable Y5LB_(0,100);
	Variable Y5UB_(0,100);





	Solver s;

	//Default Graph:

	make_Domain_Consistency(&s, &X1LB_,&X1UB_);
	make_Domain_Consistency(&s, &X2LB_,&X2UB_);
	make_Domain_Consistency(&s, &X3LB_,&X3UB_);
	make_Domain_Consistency(&s, &X4LB_,&X4UB_);
	make_Domain_Consistency(&s, &X5LB_,&X5UB_);


	make_Domain_Consistency(&s, &Y1LB_,&Y1UB_);
	make_Domain_Consistency(&s, &Y2LB_,&Y2UB_);
	make_Domain_Consistency(&s, &Y3LB_,&Y3UB_);
	make_Domain_Consistency(&s, &Y4LB_,&Y4UB_);
	make_Domain_Consistency(&s, &Y5LB_,&Y5UB_);




	make_precedence(&s,&X1LB_, &X2LB_, &X1UB_, &X2UB_, 2);
	make_precedence(&s,&X2LB_, &X3LB_, &X2UB_, &X3UB_, 2);
	make_precedence(&s,&X3LB_, &X4LB_, &X3UB_, &X4UB_, 2);
	make_precedence(&s,&X4LB_, &X5LB_, &X4UB_, &X5UB_, 2);

	make_precedence(&s,&Y1LB_, &Y2LB_, &Y1UB_, &Y2UB_, 2);
	make_precedence(&s,&Y2LB_, &Y3LB_, &Y2UB_, &Y3UB_, 2);
	make_precedence(&s,&Y3LB_, &Y4LB_, &Y3UB_, &Y4UB_, 2);
	make_precedence(&s,&Y4LB_, &Y5LB_, &Y4UB_, &Y5UB_, 2);

	//Decided disjunction :
	make_precedence(&s,&Y3LB_, &X4LB_, &Y3UB_, &X4UB_, 2);

	//Disjunctions :
	make_disjunction(&s, &X1LB_,&Y1LB_, &X1UB_, &Y1UB_,2, 2);
	make_disjunction(&s, &X2LB_,&Y2LB_, &X1UB_, &Y1UB_,2, 2);
	make_disjunction(&s, &X1LB_,&Y1LB_, &X1UB_, &Y1UB_,2, 2);
	make_disjunction(&s, &X1LB_,&Y1LB_, &X1UB_, &Y1UB_,2, 2);




	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();

	//	cout << " X  result: " << result << endl;
	//std:: cout << s << std::endl;

	//	std:: cout << X.get_solution_int_value() << std::endl;
	//	std:: cout << Y.get_solution_int_value() << std::endl;
	//	std:: cout << Z.get_solution_int_value() << std::endl;

	if(result)
	{


		cout << "\n  X :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << X1LB_.get_solution_int_value() << std::endl;
		std:: cout << X1UB_.get_solution_int_value() << std::endl;

		std:: cout << X2LB_.get_solution_int_value() << std::endl;
		std:: cout << X2UB_.get_solution_int_value() << std::endl;

		std:: cout << X3LB_.get_solution_int_value() << std::endl;
		std:: cout << X3UB_.get_solution_int_value() << std::endl;

		std:: cout << X4LB_.get_solution_int_value() << std::endl;
		std:: cout << X4UB_.get_solution_int_value() << std::endl;

		std:: cout << X5LB_.get_solution_int_value() << std::endl;
		std:: cout << X5UB_.get_solution_int_value() << std::endl;



		cout << "\n \n  Y :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << Y1LB_.get_solution_int_value() << std::endl;
		std:: cout << Y1UB_.get_solution_int_value() << std::endl;

		std:: cout << Y2LB_.get_solution_int_value() << std::endl;
		std:: cout << Y2UB_.get_solution_int_value() << std::endl;

		std:: cout << Y3LB_.get_solution_int_value() << std::endl;
		std:: cout << Y3UB_.get_solution_int_value() << std::endl;

		std:: cout << Y4LB_.get_solution_int_value() << std::endl;
		std:: cout << Y4UB_.get_solution_int_value() << std::endl;

		std:: cout << Y5LB_.get_solution_int_value() << std::endl;
		std:: cout << Y5UB_.get_solution_int_value() << std::endl;






	}
	else
		cout << "\n  UNSAT :  "<< endl;










}





void makeACounterExampleLearningwithClausesandrevisitingAdisjunction3JOBS()
{


	Variable X1LB_(0,100);
	Variable X1UB_(0,100);

	Variable X2LB_(0,100);
	Variable X2UB_(0,100);

	Variable X3LB_(0,100);
	Variable X3UB_(0,100);

	Variable X4LB_(0,100);
	Variable X4UB_(0,100);

	Variable X5LB_(0,100);
	Variable X5UB_(0,100);


	Variable Y1LB_(0,100);
	Variable Y1UB_(0,100);

	Variable Y2LB_(0,100);
	Variable Y2UB_(0,100);

	Variable Y3LB_(0,100);
	Variable Y3UB_(0,100);

	Variable Y4LB_(0,100);
	Variable Y4UB_(0,100);

	Variable Y5LB_(0,100);
	Variable Y5UB_(0,100);


	Variable Z1LB_(0,100);
	Variable Z1UB_(0,100);

	Variable Z2LB_(0,100);
	Variable Z2UB_(0,100);

	Variable Z3LB_(0,100);
	Variable Z3UB_(0,100);

	Variable Z4LB_(0,100);
	Variable Z4UB_(0,100);

	Variable Z5LB_(0,100);
	Variable Z5UB_(0,100);



	Variable KX4(1,10);
	Variable KZ3(1,10);
	Variable KZ4(1,10);


	Solver s;

	//Default Graph:

	make_Domain_Consistency(&s, &X1LB_,&X1UB_);
	make_Domain_Consistency(&s, &X2LB_,&X2UB_);
	make_Domain_Consistency(&s, &X3LB_,&X3UB_);
	make_Domain_Consistency(&s, &X4LB_,&X4UB_);
	make_Domain_Consistency(&s, &X5LB_,&X5UB_);


	make_Domain_Consistency(&s, &Y1LB_,&Y1UB_);
	make_Domain_Consistency(&s, &Y2LB_,&Y2UB_);
	make_Domain_Consistency(&s, &Y3LB_,&Y3UB_);
	make_Domain_Consistency(&s, &Y4LB_,&Y4UB_);
	make_Domain_Consistency(&s, &Y5LB_,&Y5UB_);


	make_Domain_Consistency(&s, &Z1LB_,&Z1UB_);
	make_Domain_Consistency(&s, &Z2LB_,&Z2UB_);
	make_Domain_Consistency(&s, &Z3LB_,&Z3UB_);
	make_Domain_Consistency(&s, &Z4LB_,&Z4UB_);
	make_Domain_Consistency(&s, &Z5LB_,&Z5UB_);


	make_precedence(&s,&X1LB_, &X2LB_, &X1UB_, &X2UB_, 2);
	make_precedence(&s,&X2LB_, &X3LB_, &X2UB_, &X3UB_, 2);
	make_precedence(&s,&X3LB_, &X4LB_, &X3UB_, &X4UB_, 2);
	make_precedence(&s,&X4LB_, &X5LB_, &X4UB_, &X5UB_, &KX4);

	make_precedence(&s,&Y2LB_, &Y1LB_, &Y2UB_, &Y1UB_, 2);
	make_precedence(&s,&Y1LB_, &Y3LB_, &Y1UB_, &Y3UB_, 2);
	make_precedence(&s,&Y3LB_, &Y4LB_, &Y3UB_, &Y4UB_, 2);
	make_precedence(&s,&Y4LB_, &Y5LB_, &Y4UB_, &Y5UB_, 2);

	make_precedence(&s,&Z1LB_, &Z2LB_, &Z1UB_, &Z2UB_, 2);
	make_precedence(&s,&Z2LB_, &Z3LB_, &Z2UB_, &Z3UB_, 2);
	make_precedence(&s,&Z3LB_, &Z4LB_, &Z3UB_, &Z4UB_, &KZ3);
	make_precedence(&s,&Z4LB_, &Z5LB_, &Z4UB_, &Z5UB_, &KZ4);


	//Decided disjunction :
	//make_precedence(&s,&Y3LB_, &X4LB_, &Y3UB_, &X4UB_, 2);

	//Disjunctions :
	make_disjunction(&s, &X1LB_,&Y1LB_, &X1UB_, &Y1UB_,2, 2);
	make_disjunction(&s, &X1LB_,&Z1LB_, &X1UB_, &Z1UB_,2, 2);
	make_disjunction(&s, &Z1LB_,&Y1LB_, &Z1UB_, &Z1UB_,2, 2);

	make_disjunction(&s, &X2LB_,&Y2LB_, &X2UB_, &Y2UB_,2, 2);
	make_disjunction(&s, &X2LB_,&Z2LB_, &X2UB_, &Z2UB_,2, 2);
	make_disjunction(&s, &Z2LB_,&Y2LB_, &Z2UB_, &Z2UB_,2, 2);

	make_disjunction(&s, &X3LB_,&Y3LB_, &X3UB_, &Y3UB_,2, 2);
	make_disjunction(&s, &X3LB_,&Z3LB_, &X3UB_, &Z3UB_,2, &KZ3);
	make_disjunction(&s, &Z3LB_,&Y3LB_, &Z3UB_, &Z3UB_,&KZ3, 2);

	make_disjunction(&s, &X4LB_,&Y4LB_, &X4UB_, &Y4UB_,&KX4, 2);
	make_disjunction(&s, &X4LB_,&Z4LB_, &X4UB_, &Z4UB_,&KX4, &KZ4);
	make_disjunction(&s, &Z4LB_,&Y4LB_, &Z4UB_, &Z4UB_,&KZ4, 2);

	make_disjunction(&s, &X5LB_,&Y5LB_, &X5UB_, &Y5UB_,2, 2);
	make_disjunction(&s, &X5LB_,&Z5LB_, &X5UB_, &Z5UB_,2, 2);
	make_disjunction(&s, &Z5LB_,&Y5LB_, &Z5UB_, &Z5UB_,2, 2);




	//new lower/upperbounds :


	Variable Y1LB_NEW(0,100);
	Variable Y3LB_NEW(0,100);
	Variable Y4LB_NEW(0,100);

	Variable X4UB_new(0,100);
	Variable X3UB_new(0,100);
	Variable X2UB_new(0,100);

	Variable Z4UB_new(0,100);
	Variable Z3UB_new(0,100);

	//Make propagation :

	make_lowerbound_propagation(&s, &Y1LB_NEW, &Y1LB_);
	make_lowerbound_propagation(&s, &Y3LB_NEW, &Y3LB_);
	make_lowerbound_propagation(&s, &Y4LB_NEW, &Y4LB_);

	make_upperbound_propagation(&s, &X4UB_new , &X4UB_);
	make_upperbound_propagation(&s, &X3UB_new , &X3UB_);
	make_upperbound_propagation(&s, &X2UB_new , &X2UB_);
	make_upperbound_propagation(&s, &Z4UB_new , &Z4UB_);
	make_upperbound_propagation(&s, &Z3UB_new , &Z3UB_);


	make_Domain_Consistency(&s,&Y1LB_NEW, &Y1UB_ );
	make_Domain_Consistency(&s,&Y3LB_NEW, &Y3UB_ );
	make_Domain_Consistency(&s,&Y4LB_NEW, &Y4UB_ );

	make_Domain_Consistency(&s,&X4LB_, &X4UB_new );
	make_Domain_Consistency(&s,&X3LB_, &X3UB_new );
	make_Domain_Consistency(&s,&X2LB_, &X2UB_new );
	make_Domain_Consistency(&s,&Z4LB_, &Z4UB_new );

	//DOMAIN FAILURE ON Z3!
	//make_Domain_Consistency(&s,&Z3LB_, &Z3UB_new );

	//Execute propagation :

	make_precedence(&s, &X4LB_, &Z4LB_, &X4UB_new, &Z4UB_,&KX4);

	//make_precedence X4 X3 X2

	make_precedence(&s,&X2LB_, &X3LB_, &X2UB_new, &X3UB_new, 2);
	make_precedence(&s,&X3LB_, &X4LB_, &X3UB_new, &X4UB_new, 2);

	//make precedence on Y 1 3 4

	make_precedence(&s,&Y1LB_NEW, &Y3LB_NEW, &Y1UB_, &Y3UB_, 2);
	make_precedence(&s,&Y3LB_NEW, &Y4LB_NEW, &Y3UB_, &Y4UB_, 2);


	//Propagate the two disjunctions!x2 -y2 and x1-y1

	// X2 < Y2
	s.add ((X2UB_new - Y2LB_) < 2);
	make_precedence(&s,&X2LB_, &Y2LB_, &X2UB_new, &Y2UB_, 2);

	//X1 < Y1
	make_precedence(&s,&X1LB_, &Y1LB_NEW, &X1UB_, &Y1UB_, 2);


	//Propagate the two disjunction Z4 -Y4

	s.add ((Z4UB_ - Y4LB_NEW) < 2);
	make_precedence(&s,&Z4LB_, &Y4LB_NEW, &Z4UB_new, &Y4UB_, &KZ4);

	//make precedence Z3 < z4

	make_precedence(&s, &Z3LB_,&Z4LB_, &Z3UB_new, &Z4UB_new,&KZ3);

	s.add ((Z4UB_new -Z3UB_new) == KZ3);
	s.add (Z3UB_new < Z3LB_);

	std:: cout << s << std::endl;

	Outcome result = s.depth_first_search();

	//	cout << " X  result: " << result << endl;
	//std:: cout << s << std::endl;

	//	std:: cout << X.get_solution_int_value() << std::endl;
	//	std:: cout << Y.get_solution_int_value() << std::endl;
	//	std:: cout << Z.get_solution_int_value() << std::endl;

	if(result)
	{


		cout << "\n  X :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << X1LB_.get_solution_int_value() << std::endl;
		std:: cout << X1UB_.get_solution_int_value() << std::endl;

		std:: cout << X2LB_.get_solution_int_value() << std::endl;
		std:: cout << X2UB_.get_solution_int_value() << std::endl;

		std:: cout << X3LB_.get_solution_int_value() << std::endl;
		std:: cout << X3UB_.get_solution_int_value() << std::endl;

		std:: cout << X4LB_.get_solution_int_value() << std::endl;
		std:: cout << X4UB_.get_solution_int_value() << std::endl;

		std:: cout << X5LB_.get_solution_int_value() << std::endl;
		std:: cout << X5UB_.get_solution_int_value() << std::endl;



		cout << "\n \n  Y :  "<< endl;
		cout << endl;
		cout << endl;


		std:: cout << Y2LB_.get_solution_int_value() << std::endl;
		std:: cout << Y2UB_.get_solution_int_value() << std::endl;

		std:: cout << Y1LB_.get_solution_int_value() << std::endl;
		std:: cout << Y1UB_.get_solution_int_value() << std::endl;

		std:: cout << Y3LB_.get_solution_int_value() << std::endl;
		std:: cout << Y3UB_.get_solution_int_value() << std::endl;

		std:: cout << Y4LB_.get_solution_int_value() << std::endl;
		std:: cout << Y4UB_.get_solution_int_value() << std::endl;

		std:: cout << Y5LB_.get_solution_int_value() << std::endl;
		std:: cout << Y5UB_.get_solution_int_value() << std::endl;





		cout << "\n \n  z :  "<< endl;
		cout << endl;
		cout << endl;

		std:: cout << Z1LB_.get_solution_int_value() << std::endl;
		std:: cout << Z1UB_.get_solution_int_value() << std::endl;

		std:: cout << Z2LB_.get_solution_int_value() << std::endl;
		std:: cout << Z2UB_.get_solution_int_value() << std::endl;

		std:: cout << Z3LB_.get_solution_int_value() << std::endl;
		std:: cout << Z3UB_.get_solution_int_value() << std::endl;

		std:: cout << Z4LB_.get_solution_int_value() << std::endl;
		std:: cout << Z4UB_.get_solution_int_value() << std::endl;

		std:: cout << Z5LB_.get_solution_int_value() << std::endl;
		std:: cout << Z5UB_.get_solution_int_value() << std::endl;




		cout << "\n \n  NEW VALUES !  :  "<< endl;
		cout << endl;
		cout << endl;
		cout << endl;

		cout << " Variable Y1LB_NEW(0,100);"<< Y1LB_NEW .get_solution_int_value() <<  endl;
		cout << "Variable Y3LB_NEW(0,100);"<<Y3LB_NEW.get_solution_int_value() << endl;
		cout << "Variable Y4LB_NEW(0,100);"<< Y4LB_NEW.get_solution_int_value() <<endl;

		cout << "Variable X4UB_new(0,100);"<< X4UB_new.get_solution_int_value() <<endl;
		cout << "Variable X3UB_new(0,100);"<< X3UB_new.get_solution_int_value() <<endl;
		cout << "Variable X2UB_new(0,100);"<< X2UB_new.get_solution_int_value() <<endl;

		cout << "Variable Z4UB_new(0,100); "<<Z4UB_new .get_solution_int_value() <<endl;
		cout << "Variable Z3UB_new(0,100); "<<Z3UB_new .get_solution_int_value() <<endl;



		cout << "\n \n  NEW VALUES !THE VALUES of KX4 KZ3 and KZ4  :  "<< endl;
		cout << endl;
		cout << endl;
		cout << endl;

		cout << ""<<KX4 .get_solution_int_value() <<endl;
		cout << "Variable  "<<KZ3 .get_solution_int_value() <<endl;
		cout << "Variable  "<<KZ4 .get_solution_int_value() <<endl;



	}
	else
		cout << "\n  UNSAT :  "<< endl;

}



int main(int argc, char **argv)
{

	//	testrewrite1();
	//	testrewrite2();
	//	testrewriteNotEq();
	//	testrewriteEqANDNotEq();
	//test_encode_literal();
	//	testDisjunctive();

	//	test_pointer();
	//	test_dynamic_cast();
	// test_while2();
	//makeACounterExample();

	//test_encode_latest_literal();
	check_BitSet();
	//	test_encode_latest_literal0();

	//makeACounterExampleLearningwithClausesandrevisitingAdisjunction();
	makeACounterExampleLearningwithClausesandrevisitingAdisjunction3JOBS();


}
