import autodiff as au
import numpy as np

class Minalg:

    '''
    Minimization algorithm class that works with the 
    Autodiff package.
    '''

    def __init__(self,dtrack):
        '''Initilized dtrack with the
        head node specified'''
        self.dtrack=dtrack;
    
    def Gradescent(self,start={},steps=1000):
        ''' 
        Unconstrained Gradient descent
        Variable step using the 
        '''
        if len(start)==0:
            pass
        else:
            pass
        stepsize=0.1
        varnames = list(self.dtrack.var.keys())
        for i in range(steps):
            # Compute the first differential equation
            order = self.dtrack.Topsort();
            self.dtrack.Forward(order)
            self.dtrack.Backward(order)
            gradients=np.array([i.gradient for i in self.dtrack.var.values()]);
            #Normalize and Scale
            gradients=(gradients/np.linalg.norm(gradients))*stepsize
            # Assign the newvalues
            for c,variable in enumerate(self.dtrack.var.values()):
                self.dtrack.var[varnames[c]].value=variable.value-gradients[c]
        return  {variable.name:variable.value  for variable in self.dtrack.var.values()}


        

if __name__=="__main__":
    mgraph = au.DTrack();
    x=mgraph.set_var(name="x",value=10);
    y=mgraph.set_var(name="y",value=10);
    z=mgraph.set_var(name="z",value=30);
    head_node = (x+2)**2+(y+2)**2+(z+2)**2;
    mgraph.set_header(head_node);
    Min = Minalg(mgraph);
    minimums = Min.Gradescent(steps=20000);
    print(minimums)








