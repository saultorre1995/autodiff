# Example1.
'''
Get the gradient of a function of two variables 
and plot the Graph Generated.
'''
import sys,os
from autodiff import *
# Check the result with simpy for the gradient
check=True;

# Graph Object 
graph = DTrack();
x = graph.set_var(value=5.0,name="x");
y = graph.set_var(value=3.0,name="y");
c1 = graph.set_cte(value=2.0,name="c1");
c2 = graph.set_cte(value=4.0,name="c2");
# Define the formula that will be translated as the header operator
header = (x**2*c2)+(y**2*c1)/(x*y)
# set the header for the graph
graph.set_header(header)
# sort the operators
graph.Topsort()
graph.Forward()
graph.Backward()
# Now get the gradient of the formula
print("The gradient of x is ",graph.var["x"].gradient)
print("The gradient of y is ",graph.var["y"].gradient)
# Plot the Graph
graph.PlotGraph()

if check:
    from sympy import symbols,diff;
    x,y = symbols('x y',real=True);
    dic={x:5.0,y:3.0}
    f_xy = (x**2*4.0)+(y**2*2.0)/(x*y)
    # Checking the partial differential
    df_dx = diff( f_xy , x)
    df_dy = diff( f_xy , y)
    print("The gradient of x by sympy ",df_dx.subs(dic))
    print("The gradient of y by sympy ",df_dy.subs(dic))
