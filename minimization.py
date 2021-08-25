# Written by Saul Gonzalez Resines
# Mail:sgr40@bath.ac.uk
# Bibliography : Algorithms for optimization
import autodiff as au
import numpy as np

class Minimization:
    
    '''
    Minimization algorithm 
    for the graphs created in the 
    autodifferentiation package
    '''

    def __init__(self,dtrack,variables=None,constraints=None):
        # Graph for the Minimization (with all the variables for the minimization)
        self.dtrack = dtrack;
        if variables is None:
            pass
        # Set the starting value for the variable
        elif type(variables)==dict:
            for varname in list(varibles.keys()):
                self.dtraj.var[varname].value=variables[varname];
        # Set the constraints for the calculation of the minimum
        self.constraints = constraints;
        self.output={};

class Gradescent(Minimization):
    '''
    Gradescent: Unconstrained algorithm
    '''
    def __init__(self,dtrack,variables=None,constraints=None):
        super().__init__(dtrack,variables,constraints);

    def run(self,steps,alpha=0.1):
        varnames = self.dtrack.var.keys();
        order = self.dtrack.Topsort();
        mytime=0
        for i in range(steps):
            # Calculate the gradient
            self.dtrack.Forward(order);
            self.dtrack.Backward(order);
            # Calculate the gradient of the var
            gradients = np.array([self.dtrack.var[varname].gradient for varname in varnames])
            # Normalize and scale the gradients
            gradients = (gradients/np.linalg.norm(gradients))*alpha
            # Assign the new values for the variables
            for c,varname in enumerate(varnames):
                self.dtrack.var[varname].value=self.dtrack.var[varname].value -gradients[c]
        # The return value is dictionary of the values of the minimums of the variables
        self.output = {variable.name:variable.value  for variable in self.dtrack.var.values()}

class Adagrad(Minimization):
    '''
    Adagrad minimization algorithm
    Unconstrained algorithm
    '''
    def run(self,steps,alpha=0.01):
        # eps avoid the multiplication by zero
        # alpha is the learning rate, normally 0.01
        eps = 10e-8;
        varnames = self.dtrack.var.keys();
        order = self.dtrack.Topsort();
        # Accumulative vector
        vector = np.zeros(len(varnames))
        for i in range(steps):
            self.dtrack.Forward(order);
            self.dtrack.Backward(order);
            # Calculate the gradient of the var
            gradients = np.array([self.dtrack.var[varname].gradient for varname in varnames])
            vector = vector + np.power(gradients,2);
            # Update values
            for c,varname in enumerate(varnames):
                self.dtrack.var[varname].value=self.dtrack.var[varname].value - gradients[c]*( alpha/(eps+np.sqrt(vector[c])))
        self.output = {variable.name:variable.value  for variable in self.dtrack.var.values()}


class Nesterov(Minimization):
    '''
    Nesterov minimization momentum
    Unconstrained algorithm
    '''
    pass

class RMSprop(Minimization):
    '''
    RMSprop Unconstrained Minimization
    '''
    pass

class Adam(Minimization):
    '''
    Adam Unconstrained Minimization
    '''
    pass
 


#
# 
#class Nesterov(Minimization):
#    '''
#    Nesterov algorithm
#    Unconstrained algorithm
#    '''

        

if __name__=="__main__":
    import time
    mgraph = au.DTrack();
    x=mgraph.set_var(name="x",value=10);
    y=mgraph.set_var(name="y",value=10);
    z=mgraph.set_var(name="z",value=30);
    l=mgraph.set_var(name="lambda",value=1.0);
    head_node = (x*x) + (y*y) + (z*z)-l*(x+y+z-1.0);
    #head_node = (x**4)+(y**2)+(z**2)
    mgraph.set_header(head_node);
    descent=Gradescent(mgraph)
    a=time.time()
    descent.run(steps=100000,alpha=0.01)
    #print(time.time()-a)
    #mgraph.PlotGraph()
    print(descent.output)
    #adagrad = Adagrad(mgraph)
    #adagrad.run(steps=100000,alpha=1.0)
    #print(adagrad.output)

    #Min = Minalg(mgraph);
    #minimums = Min.Gradescent(steps=5000);
    #print(minimums)








