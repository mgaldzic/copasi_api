#include "TC_api.h"
BEGIN_C_DECLS

#include "lpsolve/lp_lib.h"

TCAPIEXPORT void run(tc_matrix input) //first row = objective, rest = contraints, first two cols = arguments
{
	lprec *lp;
	int i,j,k;
	REAL * row;
	REAL * obj;
	REAL * soln;
	long o;
	double max = 0.0;
	tc_matrix output;


	/***get the stoichiometry matrix**

	all = tc_allItems();
	N = tc_getStoichiometry(all);

	if (N.rows < 1 || N.cols < 1) 
	{
	tc_errorReport("stoichiometry matrix is empty");
	return;
	}*/

	//tc_printMatrix(input);

	if (input.rows < 2 || input.cols < 3)
	{
		tc_errorReport("Incorrect input matrix");
		return;
	}

	/**initialize lpsolve***/

	if ((lp=make_lp(0,input.cols - 2)) == NULL)
	{
		tc_errorReport("lpsolve library failed to initialize");
		return;
	}

	/**setup the objective***/


	obj = malloc((input.cols - 1) * sizeof(double));
	obj[0] = 0.0;

	for (i=1; i < (input.cols-1); ++i)
	{
		obj[i] = tc_getMatrixValue(input,0,i+1);
		set_lowbo(lp,i,0.0);
		set_upbo(lp,i,100.0);
	}

	set_obj_fn(lp, obj);

	if (tc_getMatrixValue(input,0,0) > 0)
		set_maxim(lp);
	else
		set_minim(lp);

	/**setup the constraints***/

	row = malloc((input.cols-1) * sizeof(double));
	for (i=1; i < input.rows; ++i)
	{
		row[0] = 0;
		for (j=1; j < (input.cols-1); ++j)
			row[j] = tc_getMatrixValue(input,i,j+1);

		if (tc_getMatrixValue(input,i,0) == 0)
			add_constraint(lp, row, EQ , tc_getMatrixValue(input,i,1));
		else
			if (tc_getMatrixValue(input,i,0) == 1)
				add_constraint(lp, row, LE , tc_getMatrixValue(input,i,1));
			else
				add_constraint(lp, row, GE , tc_getMatrixValue(input,i,1));
	}

	free(row);

	/**run lpsolve**/

	set_outputfile(lp,"lpsolve.out");
	k = solve(lp);
	free(obj);

	print_constraints(lp, input.cols-2);
	print_objective(lp);
	print_solution(lp, input.cols-2);
	print_lp(lp);

	if (k != 0 && k != 1)
	{
		if (k < 0) tc_errorReport("lpsolve: out of memory");
		else
			if (k == 2) tc_errorReport("lpsolve: The model is infeasible");
			else
				if (k == 3) tc_errorReport("lpsolve: the model is unbounded");
				else 
					if (k == 4) tc_errorReport("lpsolve: the model is degenerative");
					else
						if (k == 5) tc_errorReport("lpsolve: numerical failure encountered");
						else
							if (k == 6) tc_errorReport("lpsolve: the abort routine returned TRUE");
							else
								if (k == 7) tc_errorReport("lpsolve: a timeout occurred.");
								else 
									if (k == 9) tc_errorReport("lpsolve: the model could be solved by presolve.");
									else
										if (k == 10) tc_errorReport("lpsolve: the B&B routine failed");
										else 
											if (k == 11) tc_errorReport("lpsolve: the B&B was stopped.");
											else
												if (k == 12) tc_errorReport("lpsolve: a feasible B&B solution was found");
												else
													if (k == 13) tc_errorReport("lpsolve: no feasible B&B solution was found");
													else
														tc_errorReport("lpsolve error");

		if (k > 0)
			delete_lp(lp);

		return;
	}
	if (k != 0)
	{
		tc_print("lpsolve: a sub-optimal solution was found.");
	}
	
	/**get solution and display on the screen***/
	soln = malloc((input.rows + input.cols - 2)*sizeof(double));
	get_primal_solution(lp, soln);


	/**output**/

	output = tc_createMatrix(1,input.cols-2);

	for (i=0; i < (input.cols-2); ++i)
		if (max < soln[input.rows+i])
			max = soln[input.rows+i];
	if (max < 1.0)
		max = 1.0;

	tc_deselect();
	for (i=0; i < (input.cols-2); ++i)
	{
		tc_setColumnName(output,i, tc_getColumnName(input,i+2));
		o = tc_find(tc_getColumnName(input,i+2));  //find the item with the column name
		if (o)
		{
			tc_displayNumber(o,soln[input.rows+i]);
			if (soln[input.rows+i] * 16.0/max < 1.0)
				tc_setLineWidth(o,1.0,0);
			else
				tc_setLineWidth(o,soln[input.rows+i] * 16.0/max,0);
			output.values[i] = soln[input.rows+i];
		}
	}

	tc_printMatrix(output);

	tc_deleteMatrix(output);
	tc_deleteMatrix(input);
	free(soln);
	delete_lp(lp);

	tc_zoom(0.99);

	return;
}
END_C_DECLS

