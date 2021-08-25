# Example3: Three minimization algorithms have been applied 
# Minibatch minimization algorithm
# Adagrad minimization algorithm
# Newton  minimization algorithm
from autodiff import *
from minimization import *
obj=DTrack();
# Set the variables and constants
x = obj.set_var(value=6.0,name="x");
y = obj.set_var(value=5.0,name="y");
c1 = obj.set_cte(value=2.0,name="c1");
# Set the header
header = (x+2.0)**2 + sin(y*c1);
obj.set_header(header);
descent=Gradescent(obj);
# alpha = learning rate
descent.run(steps=10000,alpha=0.01);
# Get a dictionary with the outputs of the minimization
print(descent.output)


