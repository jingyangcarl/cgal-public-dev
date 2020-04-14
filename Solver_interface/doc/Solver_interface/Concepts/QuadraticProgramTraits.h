class QuadraticProgramTraits

/*!
\cgalConcept

`QuadraticProgramVariable` is a concept of a variable in
the Quadratic Programming (QP) problem.

\cgalHasModel `CGAL::Variable<FT>`
*/
template <typename FT>
class QuadraticProgramVariable
{
public:
	/// \name Types
	/// @{

	/*!

	*/
	typedef unspecified_type FT;

	/*!
	A variable can be continuous, integer, or binary
	*/
	enum Variable_type { CONTINUOUS, INTEGER, BINARY };

	/// @}

	/// \name Creation 
	/// @{

	/*!
	Constructs a variable initialized with the pointer of the solver it belongs to,
	the variable type, lower bound, upper bound, name, and index.
	*/
	QuadraticProgramVariable(QuadraticProgramTraits* solver, Variable_type type, FT lb =, FT ub, const std::string& name, int idx);

	/// \name Operations 
	/// @{

	/// Returns the variable type
	Variable_type variable_type() const;

	/// Sets/Changes the variable type
	void set_variable_type(Variable_type t);

	/*!
	Returns the name of the variable.
	*/
	const std::string& name() const;

	/*!
	Sets the name of the variable.
	*/
	void set_name(const std::string& n);

	/*!
	Returns the index of the variable.
	*/
	int  index() const;

	/*!
	Sets the index of the variable.
	*/
	void set_index(int idx);

	/// Returns the solver that owns this variable
	const QuadraticProgramTraits* solver() const;
	QuadraticProgramTraits* solver();

	/// Sets the lower bound
	void set_lower_bound(FT lb);

	/// Sets the upper bound
	void set_upper_bound(FT ub);

	/// Sets both lower and upper bounds
	void set_bounds(FT lb, FT ub);

	/// Gets the lower bound
	FT lower_bound() const;

	/// Gets the upper bound
	FT upper_bound() const;

	/// Gets both lower and upper bounds
	void get_bounds(FT& lb, FT& ub) const;

	/// Gets the infinity threshold (e.g., 1e20). 
	/// Values greater than this value are considered as infinity.
	static FT infinity();

	/// Returns the value of the variable in the current solution.
	/// \note (1) Valid only if the program was successfully solved.
	///       (2) If the variable is integer and rounded == true, then the 
	///           value will be rounded to the nearest integer.
	FT solution_value(bool rounded = false) const;

	/// Sets the solution value (should be called internally by the solver).
	void set_solution_value(FT value);

	/// @}

}; /* end QuadraticProgramVariable */


/*!

\cgalConcept

`QuadraticProgramLinearConstraint` is a concept of a linear
constraint in the Quadratic Programming (QP) problem.

\cgalHasModel `CGAL::Linear_constraint<FT>`
*/
template <typename FT>
class QuadraticProgramLinearConstraint
{
public:
	/// \name Creation 
	/// @{

	/*!
	Constructs a linear constraint, initialized with the solver it belongs to,
	the lower bound, upper bound, name, and index.
	*/
	QuadraticProgramLinearConstraint(QuadraticProgramTraits* solver, FT lb, FT ub, const std::string& name, int idx);

	/// \name Operations 
	/// @{

	/*!
	Returns the name of the constraint.
	*/
	const std::string& name() const;

	/*!
	Sets the name of the constraint.
	*/
	void set_name(const std::string& n);

	/*!
	Returns the index of the constraint.
	*/
	int  index() const;

	/*!
	Sets the index of the constraint.
	*/
	void set_index(int idx);

	/// Returns the solver that owns this constraint
	const QuadraticProgramTraits* solver() const;
	QuadraticProgramTraits* solver();

	/// Sets the lower bound
	void set_lower_bound(FT lb);

	/// Sets the upper bound
	void set_upper_bound(FT ub);

	/// Sets both lower and upper bounds
	void set_bounds(FT lb, FT ub);

	/// Gets the lower bound
	FT lower_bound() const;

	/// Gets the upper bound
	FT upper_bound() const;

	/// Gets both lower and upper bounds
	void get_bounds(FT& lb, FT& ub) const;

	/// Gets the infinity threshold (e.g., 1e20). 
	/// Values greater than this value are considered as infinity.
	static FT infinity();

	/// Sets the coefficients of the constraint. 
	void  set_coefficients(const std::unordered_map<const QuadraticProgramVariable*, FT>& coeffs);

	/// Adds a coefficient to a variable of the constraint. 
	void  add_coefficient(const QuadraticProgramVariable* var, FT coeff);

	/// Returns the coefficients of the constraint. 
	const std::unordered_map<const QuadraticrProgramVariable*, FT>& coefficients() const;

	/// Gets the coefficient of the variable in this constraint. 
	FT get_coefficient(const QuadraticProgramVariable* var) const;

	/// Sets the constant term.
	void set_offset(FT value);

	/// Gets the constant term.
	FT offset() const;

	/// Clears all variables and sets the constant term to zero.
	/// Useful to reuse the object to define a new linear constraint.
	void clear();

	/// @}

}; /* end QuadraticProgramLinearConstraint */


/*!

\cgalConcept

`QuadraticProgramLinearObjective` is a concept of the linear
objective function in the Quadratic Programming (QP) problem.

\cgalHasModel `CGAL::Linear_objective<FT>`
*/
template <typename FT>
class QuadraticProgramLinearObjective
{
public:
	/// \name Types
	/// @{

	/// The objective sense (i.e., optimization direction)
	enum Sense { MINIMIZE, MAXIMIZE, UNDEFINED };

	/// @}

	/// \name Creation 
	/// @{

	/*!
	Constructs a linear objective, initialized with the solver it belongs to
	and the objective sense.
	*/
	QuadraticProgramLinearObjective(QuadraticProgramTraits* solver, Sense sense);

	/// \name Operations 
	/// @{

	/// Sets the objective sense.
	void  set_sense(Sense sense);

	/// Gets the objective sense.
	Sense sense() const;

	/// Sets the coefficients of the constraint. 
	void  set_coefficients(const std::unordered_map<const QuadraticProgramVariable*, FT>& coeffs);

	/// Adds a coefficient to a variable of the constraint. 
	void  add_coefficient(const QuadraticProgramVariable* var, FT coeff);

	/// Returns the coefficients of the constraint. 
	const std::unordered_map<const QuadraticProgramVariable*, FT>& coefficients() const;

	/// Gets the coefficient of the variable in this constraint. 
	FT get_coefficient(const QuadraticProgramVariable* var) const;

	/// Sets the constant term.
	void set_offset(FT value);

	/// Gets the constant term.
	FT offset() const;

	/// Clears the objective (i.e., removes all variables, resets the 
	/// objective sense to UNDEFINED). Useful to reuse the object 
	/// to define a new linear objective.
	void clear();

	/// @}

}; /* end QuadraticProgramLinearObjective */


/*!
\ingroup PkgSolverConcepts
\cgalConcept

@brief Concept describing the set of requirements for (constrained or unconstrained)
Quadratic Programming (QP) problems. A model of this concept stores the integer
variables, linear objective, and linear constraints (if any) and provides a method
to solve the problem.

\cgalHasModel `CGAL::Shape_regularization::CGAL_solver`
\cgalHasModel `CGAL::Shape_regularization::OSQP_solver`
*/
template <typename FT>
class QuadraticProgramTraits
{
public:

	/// \name Creation
	/// @{

	/*!
	Default constructor.
	*/
	QuadraticProgramTraits();

	/// @}

	/// \name Operations
	/// @{

	/// Creates a single variable, adds it to the solver, and returns its pointer.
	/// \note Memory is managed by the solver and will be automatically released.
	QuadraticProgramVariable* create_variable(Variable_type type, FT lb, FT ub, const std::string& name);

	/// Creates a set of variables, adds them to the solver, and returns their pointers.
	/// \note (1) Variables will be given default names, e.g., x0, x1...; 	
	///		  (2) Memory is managed by the solver and will be automatically released when the
	///		  solver is destroyed.
	std::vector<QuadraticProgramVariable*> create_variables(std::size_t n);

	/// Creates a single linear constraint, adds it to the solver, and returns the pointer.
	/// \note Memory is managed by the solver and will be automatically released when the
	///		  solver is destroyed.
	QuadraticProgramLinearConstraint* create_constraint(FT lb, FT ub, const std::string& name);

	/// Creates a set of linear constraints, adds them to the solver, and returns their pointers.	
	/// \note (1) Constraints will be given default names, e.g., c0, c1...
	///		  (2) Memory is managed by the solver and will be automatically released when the
	///		  solver is destroyed.
	std::vector<QuadraticProgramLinearConstraint*> create_constraints(std::size_t n);

	/// Creates the objective function and returns the pointer.
	/// \note Memory is managed by the solver and will be automatically released when the
	///		  solver is destroyed.
	QuadraticProgramLinearObjective* create_objective(Sense sense);

	/// Returns the number of variables
	std::size_t number_of_variables() const;

	/// Returns the variables
	const std::vector<QuadraticProgramVariable*>& variables() const;
	std::vector<QuadraticProgramVariable*>& variables();

	/// Returns the number of constraints
	std::size_t number_of_constraints() const;

	/// Returns the constraints
	const std::vector<QuadraticProgramLinearConstraint*>& constraints() const;
	std::vector<QuadraticProgramLinearConstraint*>& constraints();

	/// Returns the objective
	const QuadraticProgramLinearObjective * objective() const;
	QuadraticProgramLinearObjective * objective();

	/// Solves the program. Returns false if failed.
	bool solve();

	/// Returns the result. 
	/// \note (1) Result is valid only if the solver succeeded.
	///       (2) Each entry in the result corresponds to the variable with the
	///			 same index in the program.
	const std::vector<FT>& solution() const;

	/// Returns the error message.
	/// \note This function should be called after call to solve().
	const std::string& error_message() const { return error_message_; }

	/// Clears all variables, constraints, and the objective.
	void clear();

	/// @}

}; /* end QuadraticProgramTraits */
