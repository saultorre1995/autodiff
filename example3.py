# Example3: Three minimization algorithms have been applied 
'''
First Order Minimization Algorithms Implementation:

Gradient Descent minimization algorithm
Adagrad minimization algorithm
Nesterov minimization algorithm
Adadelta minimization algorithm
Adam Minimizartion algorithm
'''

from autodiff import *
from minimization import *

def set_vars(graph,dic):
    for element in dic:
        graph.var[element].value=dic[element];


obj=DTrack();
dic_start={"x":6.0,"y":5.0}
# Set the variables and constants
x = obj.set_var(value=6.0,name="x");
y = obj.set_var(value=5.0,name="y");
c1 = obj.set_cte(value=2.0,name="c1");

# Set the header
header = (x+2.0)**2 + sin(y*c1);
obj.set_header(header);

# Gradient Descent
descent=Gradescent(obj);
# alpha = learning rate
descent.run(steps=5000,alpha=0.01);
# Get a dictionary with the outputs of the minimization
print("Gradient Descent Output",descent.output)

# Adagrad 
set_vars(obj, dic = dic_start)
adagrad = Adagrad(obj)
adagrad.run(steps=10000,alpha=0.1);
print("Adagrad Output",adagrad.output)

# Nesterov 
set_vars(obj, dic = dic_start)
nesterov = Nesterov(obj)
# Beta Momentum decay
nesterov.run(steps=5000,alpha=0.01,beta=0.1)
print("Nesterov Output",nesterov.output)

# Adadelta
set_vars(obj, dic = dic_start);
adadelta = Adadelta(obj);
adadelta.run(1000,beta_grad=0.5,beta_diff=0.2,eps=0.0001)
print("Adadelta Output",adadelta.output)

# Adam
set_vars(obj ,dic = dic_start);
adam = None





