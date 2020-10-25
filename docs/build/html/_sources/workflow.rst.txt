Example workflow with LoopDetect
################################


Installation
------------

LoopDetect is on PyPI and can be installed via pip in the command shell (you might have to use pip3 if using Python3).

.. code-block:: bash

	pip install loopdetect

Within Python, import the Loopdetect's core module like this in order to use LoopDetect detection functions:

.. code-block::

	import loopdetect.core 

If you want to use example functions, import the example module like this:

.. code-block::

	import loopdetect.examples 



In brief and quick start
------------------------

The package LoopDetect enables determining all feedback loops of
an ordinary differential equation (ODE) system at user-defined values of 
the model parameters and of the modelled variables.

The following calls within Python report (up to 10) feedback loops for an ODE system 
determined by a function (here the example function `func_POSm4` provided along with LoopDetect) at variable values `s_star` (here, these are all equal to 1). Additional arguments to the 
example function are supplied. 

.. code-block:: python

	# Import LoopDetect core and numpy
	import numpy as np
	import loopdetect.core as ld
	import loopdetect.examples as lde
	# An example ODE system with function func_POSm4 is provided in 
	# loopdetect.examples, it has 4 variables.
	# Variable values for the 4 variables (as tuple, cast into a list)
	s_star = [(1,1,1,1)]
	# Define further arguments of func_POSm4
	klin = np.ones(1,8)
	knonlin = (2.5,3)
	# compute loops
	res_tab = ld.find_loops_vset(lde.func_POSm4,vset=s_star,klin=klin,
                           knonlin=knonlin,max_num_loops=10)
	# The loop list, a pandas dataframe, is accessed like this:
	res_tab['loop_rep'][0]

The following describes how information in the output pandas dataframe
relates to information on selected loops.

.. code-block:: python

	# Retrieve only the second loop of the list. 
	res_tab['loop_rep'][0].iloc[1]
	# It is a positive feedback loop (sign in the loop list equals +1) 
	# of length 3 in that variable 2 (index 1) regulates variable 3 
	# (index 2), variable 3 regulates variable 4 (index 3), and 
	# variable 4 regulates variable 2. 
	# This returns the corresponding signed Jacobian matrix for 
	# the whole system.
	res_tab['jac_rep'][0]


Introduction
------------

Ordinary differential equation (ODE) models are used frequently to
mathematically represent biological systems. Feedback loops are important 
regulatory features of biological systems and can give rise to different
dynamic behavior such as multistability for positive feedback loops or
oscillations for negative feedback loops.

The feedback loops in an ODE system can be detected with the help of its 
Jacobian matrix, the matrix of partial derivatives of the variables. 
It captures all interactions between the variables and gives rise to the
interaction graph of the ODE model. In this graph, each modelled variable
is a node and non-zero entries in the Jacobian matrix are (weighted) 
edges of the graph. Interactions can be positive or negative, according 
to the sign of the Jacobian matrix entry.

Directed path detection in this graph is used to determine all feedback 
loops (in graphs also called cycles or circuits) of the system. They are
marked by a set of directed interactions forming a chain in which only 
the first and the last node (variable) is the same. Thereby, self-loops 
(loops of length one) can also occur.

LoopDetect allows for detection of all loops of the graph and also reports
the sign of each loop, i.e. whether it is a positive feedback loop
(the number of negative interactions is even) or a negative feedback loop 
(the number of negative interactions is uneven).
The output is a table that captures the order of the variables forming 
the loop, their length and the sign of each loop. 

Jacobian determination in LoopDetect relies on the Python module `numdifftools`, and path 
finding in graphs uses algorithms supplied in the Python module `networkx`. The output loop lists rely on the pandas dataframe format.


Solving the ODE model to generate variable values of interest
-------------------------------------------------------------

Solving an ODE model can be performed with the `odeint` function from `scipy.integrate`.
Note: You can skip this step if you already have a point of interest in state 
space, or if you want to use dummy values for the variables such as 
`s_star = (1,2,3,4)`.

.. code-block:: python
	
	# Import required modules 
	from scipy.integrate import odeint 
	import loopdetect.examples as lde
	# We use the example ODE system with function func_POSm4, positive 
	# feedback chain model from [Baum et al., 2016], 4 variables.
	# It is supplied as function in the examples module in LoopDetect.
	# Kinetic parameters of the model, further arguments to func_POSm4
	klin = (165,0.044,0.27,550,5000,78,4.4,5.1)
	knonlin = (0.3,2)
	# Define a helper function: odeint requires a function depending 
	# only on the variables (first argument) and the time, t, but 
	# func_POSm4 is independent from t and carries two more arguments 
	# (klin, knonlin).
	def func_POSm4_help(x,t):
		return(lde.func_POSm4(x,klin,knonlin))
	# Solve the system using odeint at time points 0, 1, ..., 100, 
	# initial vector: (1,2,3,4)
	sol = odeint(func_POSm4_help, y0 = (1,2,3,4), 
		t = np.linspace(0,100,11))
	# We set the last point of the numeric solution as point of interest
	s_star = sol[-1]


State variable values of interest could be steady state values,
values at a specific point in time (e.g. after a stimulus) or even multiple sets 
of values (see section *Determining loops over multiple sets of variable values*). 


Calculating the Jacobian matrix
-------------------------------

The function `Jacobian` from the numdifftools package can be used to determine 
numerically the Jacobian matrix of an ODE system at a certain set of values 
for the variables, `s_star`. 
The approach is that of finite differences (with real step) or 
complex step approach, the latter of which is supposed to deliver more exact 
results. Please note that if the complex-step approach should be performed, the supplied function for that the Jacobian is to be determined has to deliver complex values and is only allowed to employ differentiable functions. Please check formats in the section *Example code for a possible input function*.

The input function, in the example `func_POSm4_comp` (positive feedback chain model 
from [Baum et al., 2016] with complex values), defines the time derivatives of the modelled variables as a vector:  :math:`f_i(S)=dS_i/dt`. Note that only those input arguments to the 
function that encode the modelled variables (and hence in whose direction the 
partial derivatives are taken) are allowed to be called `x`, they should always go as first argument to the function to derive.

.. code-block:: python
	
	# Import the required packages. 
	import numpy as np
	import numdifftools
	import loopdetect.examples as lde
	# Define kinetic parameters of the model that are arguments of the 
	# example function func_POSm4_comp provided in LoopDetect
	klin = tuple(165,0.044,0.27,550,5000,78,4.4,5.1)
	knonlin = tuple(0.3,2)
	# Note that we defined s_star in the section above. It could also 
	# be defined as simple tuple, s_star=(1,2,3,4).
	j_matrix = numdifftools.Jacobian(lde.func_POSm4_comp,method="complex")(s_star,
                               klin=klin,knonlin=knonlin).real                          
	signed_jacobian = np.sign(j_matrix)


The (i,j)th entry of the Jacobian matrix denotes the partial derivative
of variable :math:`S_i` with respect to variable :math:`S_j`, 
:math:`J_{ij}=\delta S_i/\delta S_j`, which is positive if :math:`S_j` has a direct positive 
effect on :math:`S_i`, negative if :math:`S_j` has a direct negative effect on :math:`S_i` and zero
if :math:`S_j` does not have a direct effect on :math:`S_i`. For example, the entry in row 2,
column 4, :math:`J_{24}`, of the `signed_jacobian` matrix above (`signed_jacobian[1,3]`) is positive, meaning that
in the underlying ODE model, variable 4 positively regulates variable 2.


Computing all feedback loops
----------------------------

The Jacobian matrix is used to compute feedback loops in the generated 
interaction graph. The default function for this is `find_loops`, in that 
strongly 
connected components are determined to reduce runtime. For smaller
systems, the function `find_loops_noscc` skips this step and thus can be
faster. The optional second input argument, `max_num_loops`, sets an upper 
limit to the number of detected and reported loops and thus can prevent overly 
long runtime (but also potentially not all loops are returned).

.. code-block:: python

	# Import the required package
	import loopdetect.core as ld
	# Determine the loop_list from Jacobian matrix j_matrix
	loop_list = ld.find_loops(j_matrix)
	# The signed Jacobian matrix can be supplied instead, 
	# delivering the same results.
	loop_list = ld.find_loops(signed_jacobian)


The output loop list is a pandas dataframe with one row for each detected loop. 
The column 

	- 'loop' contains the order of the variables that form the loop (as tuple); 
	- 'length' contains the loop length (i.e. the number of variables involved);
	- 'sign' denotes whether the loop is negative, -1, or positive, 1.

Finding certain loops or edges in a loop list
---------------------------------------------

The dataframe can be queried as usual for single loops or loops of a certain 
length or sign. 

.. code-block:: python 
	
	# Retrieve second to fourth loop and its properties
	loop_list.iloc[1:4,] 
	# In order to obtain only the vector of the order of species of a 
	# loop, you have to call the first, the 'loop', column.
	loop_list.iloc[5,0]
	# Retrieve all loops of length 4
	loop_list.loc[loop_list['length']==4] 

The LoopDetect function `loop_summary` provides a convenient report on total 
number of loops, subdivided by their lengths (columns) and signs (`total`, `pos`,
`neg`) in a pandas dataframe.

.. code-block:: python 

	lsum = loop_summary(loop_list)
	#retrieve the loop counts for all loops of length 3
	lsum[3]
	#loop counts subdivided by length for all negative feedback loops
	lsum.loc['neg']
	#retrieve all negative loops of length 1
	lsum[1]['neg']

One can filter the loop list for loops containing specific variables, for 
example the one with index 2:

.. code-block:: python
	
	#define node (variable) of interest
	noi = 2
	#only return loops containung this variable
	loop_list[[noi in loop for loop in loop_list['loop']]]

The LoopDetect function `find_edge` can be used to search a loop list for loops 
containing specific edges defined by the indices of the ingoing and outgoing 
nodes. This example returns the indices of all loops with a regulation of node 
with index 3 by node with index 2. These are only two here.

.. code-block:: python

	# Obtain the indices of the loops with edge '2 regulates 3'
	loop_edge_ind = find_edge(loop_list,source_node=2,target_node=3)
	loops_with_edge_2_to_3 = loop_list[loop_edge_ind]

Writing to files and loading loop lists
---------------------------------------

For input/output procedures, all Python pandas functionalities can be used.
Loop lists can be saved and loaded as Python objects using the pandas pickling functions. 

.. code-block:: python

	# this writes a Python-specific pickled file 'looplist_func_POSm4'
	loop_list.to_pickle('looplist_func_POSm4')
	# this loads the loop list back from the file. We have to import 
	# pandas in order to perform this.
	import pandas as pd
	ll1 = pd.read_pickle('looplist_func_POSm4')

Saving the loop list into tab-separated table files and loading them back is also possible with pandas table writing options.

.. code-block:: python

	# this saves the loop list as tab-separated file (without index)
	loop_list.to_csv('looplist_func_POSm4.tsv',sep="\t",index=False)
	# this reads the tsv back as loop list (with new indices)
	# pandas has to be imported (see above)
	ll2 = pd.read_csv('looplist_func_POSm4.tsv',sep="\t")


Computing feedback loops over multiple sets of variable values of interest
--------------------------------------------------------------------------

In this example of a model of the bacterial cell cycle [Li et al., 2008],
it is demonstrated how feedback loops can be determined over multiple sets of
variable values. Here, it is focused on the solution of the ODE
systems along the time axis (provided as data in the package). 

For convenience, the example ODE model function func_li08 is provided as a 
LoopDetect function,
the kinetic parameters are defined within the function. The function returns 
a vector of the time derivatives (in the same order as the modelled variables 
in the arguments).

.. code-block:: python
	
	# import required packages
	import loopdetect.core as ld
	import loopdetect.examples as lde
	# load solution over time as set of variable values of interest
	sols = lde.load_li08_sol()
	# convert the solution to a list of tuples, remove the time column
	sols_as_tuples = [tuple(sols.iloc[i,1:20]) for i in range(len(sols))]
	# compute all different loop lists along the solution
	# attention - this may take several minutes due to computing the 
	# Jacobians at each set of variable values
	loop_results = ld.find_loops_vset(lde.func_li08,vset=sols_as_tuples,
		numdiff_method='central',max_num_loops=100000,
		compute_full_list=False,t=0)

The solutions of the example ODE model give rise to seven different loop lists 
that are saved as elements of the list `loop_results['loop_rep']`.
Here, an example of a resulting loop list is given.

.. code-block:: python
	
	# this is the 6th loop list
	loop_list_1 = loop_results['loop_rep'][5]
	# here we show only loops of length > 1 (i.e. no self-loops)
	loop_list_1.loc[loop_list_1['length']>1]

The entry `loop_results['loop_rep_index']` is a vector of the same length as the 
number of different states at which the loops were determined 
(`len(sols_as_tuples)`), and returns which entry of `loop_results['loop_rep']` belongs
to each input state.

.. code-block:: python

	# determine which sets of variables belong to the 6th loop list 
	# from above; yields a vector of the same size as sols_as_tuples
	loop_results['loop_rep_index']==5
	# the indices can be retrieved using the numpy function `where`
	np.where(loop_results['loop_rep_index']==5)

Similarly, `loop_results['jac_rep']` and `loop_results['jac_rep_index']` capture the different
Jacobian matrices and for each input state which Jacobian belongs to it.

Results of this analysis could be plotted along the solution and analyzed
further to discover reasons of changing loops. Please note that in order 
to obtain the sample solution that can be loaded with `ld.load_li08_sol()`, 
also event functions
are required; the solution cannot be retrieved from integrating 
`func_li08` alone. Please refer to the model's publication [Li et al.,
2008] for details.


Comparing two loop lists
------------------------

LoopDetect provides a function for comparing the loops of two systems, 
`compare_loop_list`. For a meaningful comparison, the loop indices in
the compared systems should point to the same variables. This could be the case 
when regulations change within one system between different sets of variables of 
interest (along a dynamic trajectory, at different steady states of the system),
or when comparing different systems in which one or more regulations are altered 
(as in the below example the positive feedback chain model vs. the negative 
feedback chain model, [Baum et al., 2016]).

We first compute the list of feedback loops for the first model (positive feedback chain model).

.. code-block:: python
	
	# import LoopDetect and numdifftools
	import loopdetect.core as ld
	import loopdetect.examples as lde
	import numdifftools 
	# Set kinetic parameters
	klin = (165,0.044,0.27,550,5000,78,4.4,5.1)
	knonlin = (0.3,2)
	# Compute the Jacobian matrix of the system at the state (1,1,1,1)
	j_matrix = numdifftools.Jacobian(lde.func_POSm4_comp,
		method="complex")((1,1,1,1),klin,knonlin).real
	# Compute all loops for this Jacobian
	loop_list_pos = ld.find_loops(j_matrix)

Second, we compute the list of feedback loops for a second model that is slightly altered with respect to the first model; in particular, it carries the same variables (negative feedback chain model).

.. code-block:: python

	# The altered regulation affects two entries of the Jacobian matrix. 
	# Parameter values and the set of variable values remain identical.
	import numpy as np
	j_matrix_neg = np.sign(j_matrix)
	j_matrix_neg[0:2,3] = -np.sign(j_matrix[0:2,3])
	# Compute the loop list for this Jacobian with altered regulation
	loop_list_neg = ld.find_loops(j_matrix_neg);

Now, we compare the loop lists with the dedicated LoopDetect function `compare_loop_list`. It returns a dictionary with different numpy arrays as entries that are described below.

.. code-block:: python
	
	# Compute comparison
	loop_comparison = ld.compare_loops(loop_list_pos,loop_list_neg)
	# Only the four self-loops remain identical in both systems. 
	# Their indices with respect to the first input list, loop_list_pos, 
	# are saved in ind_a_id.
	loop_list_pos.loc[loop_comparison['ind_a_id']]
	# ind_b_id saves the indices of the corresponding loops in the 
	# second loop list, loop_list_neg.
	loop_list_neg.loc[loop_comparison['ind_b_id']]
	# Two loops are the same in both systems but they have switched 
	# their signs. Their indices in the first loop list are saved in 
	# ind_a_switch.
	loop_list_pos.loc[loop_comparison['ind_a_switch']]
	# Their indices in the second loop list are saved in ind_b_switch. 
	loop_list_neg.loc[loop_comparison['ind_b_switch']]
	# All loops in the first system do also occur in the second system, 
	# i.e. the entry ind_a_notin that gives the indices of these loops
	# within the first loop list is an empty array.
	loop_comparison['ind_a_notin']



Example for a possible input function
-------------------------------------

We here show how a function that could be an input to `find_loops_vset`.
It is the function that is supplied within the example module of the LoopDetect package, `func_POSm4`, implementing 
a chain of length 4 with the last species feeding back on the conversion between the first and 
the second species.

.. literalinclude:: func_POSm4.py
	:lines: 3, 32-37


In order to use this function with complex-step derivative, we simply set the output type to a complex numpy array. This requires numpy and should be a sufficient procedure for most other functions as well.

.. code-block:: python

	dx = numpy.zeros(4,dtype='complex')

Please note that the performed operations within the function should support complex numbers
and be differentiable!




References
----------

Baum K, Politi AZ, Kofahl B, Steuer R, Wolf J. Feedback, Mass 
Conservation and Reaction Kinetics Impact the Robustness of Cellular 
Oscillations. PLoS Comput Biol. 2016;12(12):e1005298.

Li S, Brazhnik P, Sobral B, Tyson JJ. A Quantitative Study of the 
Division Cycle of Caulobacter crescentus Stalked Cells. Plos Comput Biol. 
2008;4(1):e9.



