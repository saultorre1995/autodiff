# Autodifferentiation Library
# Written by Saul Gonzalez Resines
# Mail:sgr40@bath.ac.uk

import warnings
import numpy as np

def sin(a,owner=None):
    '''
    Sine Function calling an inner class of the DObject that
    wraps the operator into an Operator object.
    The input must be a DObject or a float,integer but in the last 
    case you need to include the owner DTrack of the Operator.
    The a must be in Radians.
    '''
    if isinstance(a,DObject):
        return a.sin();
        pass;
    elif isinstance(a,float) or isinstance(a,int) and isinstance(owner,DTrack):
        constant=owner.set_cte(value=a);
        return constant.sin();
    else:
            raise TypeError("Incompatible types, or the object DTrack not defined.")

def cos(a,owner=None):
    '''
    Cosine Function calling an inner class of the DObject that
    wraps the operator into an Operator object.
    The input must be a DObject or a float,integer but in the last 
    case you need to include the owner DTrack of the Operator.
    The a must be in Radians.
    '''
    if isinstance(a,DObject):
        return a.cos();
    elif isinstance(a,float) or isinstance(a,int) and isinstance(owner,DTrack):
        constant=owner.set_cte(value=a);
        return constant.cos();
    else:
            raise TypeError("Incompatible types, or the object DTrack not defined.")

def logb (a,b,owner=None):
    '''
    Log base 10 calling an inner class of the DObject that
    wraps the operator into an Operator object.
    The input must be a DObject or a float,integer but in the last 
    case you need to include the owner DTrack of the Operator.
    '''
    if isinstance(a,DObject):
        return a.log10();
    elif isinstance(a,float) or isinstance(a,int) and isinstance(owner,DTrack):
        constant=owner.set_cte(value=a);
        return constant.log10();
    else:
            raise TypeError("Incompatible types, or the object DTrack not defined.")


def logn (a,owner=None):
    '''
    Natural Logarithm calling an inner class of the DObject that
    wraps the operator into an Operator object.
    The input must be a DObject or a float,integer but in the last 
    case you need to include the owner DTrack of the Operator.
    '''
    if isinstance(a,DObject):
        return a.logn();
    elif isinstance(a,float) or isinstance(a,int) and isinstance(owner,DTrack):
        constant=owner.set_cte(value=a);
        return constant.logn();
    else:
            raise TypeError("Incompatible types, or the object DTrack not defined.")


#def exp (a,owner=None):
#    '''
#    Exponentiation
#    '''
#    pass;
#
#def sqrt (a,owner=None):
#    '''
#    Square Root
#    '''
#    pass;


    

class DObject:
    """
    A class that encodes all the
    objects (variables,constants,operators)
    that are supported by the package.
    """
    def __init__(self,owner=None):
        self.owner=owner
        pass
    # Operator overloading over the class DObject
    # for the sum and the multiplication.

    def _wrapper(self,operator,other):
        # Just create the object
        if self.owner is None:
            print("The variable with name "+self.name+" does not have any owner")
            return None;
        else:
            if isinstance(other,DObject):
                return self.owner.set_ope(operator,self,other);
            elif isinstance(other,float) or isinstance(other,int):
                constant = self.owner.set_cte(other);
                return self.owner.set_ope(operator,self,constant)
            else:
                raise TypeError("Incompatible Types");
    
    def _wrapper_single(self,operator):
        return self.owner.set_ope(operator,self,None)
    
    def sin(self):
        return self._wrapper_single(Sine);

    def cos(self):
        return self._wrapper_single(Cosine);

    def logn(self):
        return self._wrapper_single(LogN);

    #def cos(self):
    #    return self._wrapper_single(Cosine);

    def __add__(self,other):
        return self._wrapper(Add,other);

    def __sub__(self,other):
        return self._wrapper(Sub,other);

    def __mul__(self,other):
        return self._wrapper(Multiply,other);

    def __pow__(self,other):
        return self._wrapper(Power,other);
    
    def __truediv__(self,other):
        return self._wrapper(Divide,other);
    
    #def sin(self):
    #    return self._
    

class DTrack:
    # This is the object the get all the DObjects
    def __init__(self,head_node=None):
        self.head_node=head_node;
        self.cte={};
        self.var={};
        self.ope={};
        ordering=[];
        # Numbering the different  operators ad DObjects
        self.n_cte  = 0;
        self.n_var  = 0;
        self.n_add  = 0;
        self.n_sub  = 0;
        self.n_mul  = 0;
        self.n_pow  = 0;
        self.n_sine = 0;
        self.n_cos  = 0;
        self.n_div  = 0;

    def __repr__(self):
        return "DTrack Object"

    def set_var(self,value,name=None):
        # Set Variable
        this = Variable(self,value,name);
        return this;
        
    def set_cte(self,value,name=None):
        # Set Constant
        this = Constant(self,value,name);
        return this;

       
    def set_ope(self,operator,a,b):
        # Set Operator
        this=operator(self,a,b);
        return this;

    def set_header(self,head):
        # Set the function or header
        self.head_node=head;

    def get_var(self,var):
        if isinstance(Variable):
            return self.var[var.name]
        elif isinstance(str):
            return self.var[var]
        else:
            raise TypeError("Not a Variable or string");
        return None


    def reset(self):
        # Restart the full of the instances (variables)
        self.head_node=None;
        self.cte.clear();
        self.var.clear();
        self.ope.clear();
        self.ordering=[];
        self.n_cte=0;
        self.n_var=0;
        self.n_add=0;
        self.n_sub=0;
        self.n_mul=0;
        self.n_pow=0;
        self.n_div=0;
        self.n_sine=0;
        self.n_cos=0;

        
    def Topsort(self):
        vis = set();
        order = [];
        def _check(node):
            if node not in vis:
                # avoid the repeat of the same variable or constant;
                # operators are going always to be different;
                vis.add(node);
                if isinstance(node,Operator):
                    for inp in node.inputs:
                        _check(inp);
                order.append(node);

        if self.head_node is None:
            print("The head_node is None so not topological order executed");
            return None
        else:
            _check(self.head_node)
            # set ordering
            self.ordering=order
            return order;
    
    def Forward(self,order=[]):
        if len(order)==0:
            order=self.ordering;
        for node in order:
            # Add the placeholder option for an easier handling
            if isinstance(node,Operator):
                node.value = node.forward(*[inp_node.value for inp_node in node.inputs]);

    def Backward(self,order=[]):
        if len(order)==0:
            order=self.ordering;
        # The last of the head_node gradient is always one
        order[-1].gradient=1;
        vis=set();
        for node in reversed(order):
            if isinstance(node,Operator):
                inputs = node.inputs;
                grads  = node.backward(*[inp.value for inp in inputs],dout=node.gradient)
                for inp,grad in zip(inputs,grads):
                    if inp not in vis: 
                        inp.gradient = grad;
                    else:
                        inp.gradient += grad;
                    vis.add(inp)
        return [node.gradient for node in order]


    def Forward_Obj(self,order=[]):
        if len(order)==0:
            order=self.ordering;
        for node in order:
            if isinstance(node,Operator):
                node.value_obj = node.forward(*node.inputs)
    
    def Backward_Obj(self,order=[]):
        if len(order)==0:
            order=self.ordering;
        # The last of the head_node gradient is always one
        vis=set();
        order[-1].gradient_obj = self.set_cte(value=1.0);
        for node in reversed(order):
            if isinstance(node,Operator):
                grads  = node.backward(*node.inputs,dout=node.gradient_obj)
                for inp,grad in zip(node.inputs,grads):
                    if inp not in vis: 
                        inp.gradient_obj = grad;
                    else:
                        inp.gradient_obj += grad;
                    vis.add(inp)
        return [node.gradient_obj for node in order]

    def GetGradObj(self):
        self.Forward_Obj();
        self.Backward_Obj();
        return None
    
    def Simplify(self):
        '''Try to reduce 
        the ammount of nodes 
        for the fast computation.
        Eliminate the nodes of two constants.
        '''
        pass
        
    
    def PlotGraph(self):
        ''' Plots the graph 
        starting for the headnode
        '''
        import matplotlib as mpl
        import matplotlib.pyplot as plt 
        if self.head_node is None:
            raise TypeError("The head_node is equal to None. Maybe it was not set.")
            return ;
        # Create The Figure and a single axis
        fig = plt.figure();
        ax = fig.add_subplot(1,1,1);
        start_xy=(0,0)
        polygons=[];
        maxnodes=[];
        def _checknode(node,count=0):
            if isinstance(node,Operator):
                count+=1;
                for inp in node.inputs:
                    _checknode(inp,count);
            else:
                maxnodes.append(count)
        
        # Here execute the code for checking the maximum distance
        _checknode(self.head_node,0);
        totlev=max(maxnodes);
        maxdist=3**(max(maxnodes));

        def _plotnode(node,pos,level=1):
            # Function to plot the operators with a circle
            size=10/level

            dist_y = maxdist/(level);
            dist_x = 6*maxdist/totlev;
            # Get the info
            if isinstance(node,DObject):
                if node.value is None:
                    valnode=""
                else:
                    valnode="{:.2f}".format(node.value)
                if node.gradient is None:
                    gradval=""
                else:
                    gradval="{:.2f}".format(node.gradient)
                info = node.name+"\n"+"val: "+valnode+\
                        "\n"+"grad: "+\
                        gradval
            # Plot the info            
            if isinstance(node,Operator):
                polygons.append(mpl.patches.Circle(pos,radius=1,color="lime"))
                ax.annotate(info,xy=pos,ha="center",va="center",fontsize=size);
                level+=1;
                if len(node.inputs)==2:
                    # if the operator has two inputs
                    post=(pos[0]+dist_x, pos[1]+dist_y);
                    polygons.append(mpl.patches.Polygon(([pos,post])))
                    _plotnode(node.inputs[0],post,level);
                    post=(pos[0]+dist_x, pos[1]-dist_y);
                    polygons.append(mpl.patches.Polygon([pos,post]))
                    _plotnode(node.inputs[1],post,level);
                elif len(node.inputs)==1:
                    # if the operator has another input
                    post=(pos[0]+dist_x,pos[1]+dist_y); 
                    _plotnode(node.inputs[0],post,level);
                else:
                    # Do not include as not possible
                    warnings.warn("The input of the operator "+node.name+" is "+\
                            str(len(node.inputs))+" which is a case not implemented in PlotGraph."+\
                            "Some nodes will be missing.")
                
            elif isinstance(node,Variable):
                polygons.append(mpl.patches.Circle(pos,radius=1,color="red"))
                ax.annotate(info,xy=pos,ha="center",va="center",fontsize=size); 
            else:
                polygons.append(mpl.patches.Circle(pos,radius=1,color="cyan"))
                ax.annotate(info,xy=pos,ha="center",va="center",fontsize=size);
        
        # Maximum distance required for the tree
        _plotnode(self.head_node,start_xy,level=1);
        # Patch the lines and the circles to the plot
        for pol in polygons:
            ax.add_patch(pol)
        ax.axis("equal")
        ax.autoscale_view();
        plt.tight_layout()
        plt.show()
        # Return the Figure Object
  

    def GetGrad(self,which=[]):
        if len(self.ordering)==0:
            self.Topsort();
        self.Forward();
        grads=self.Backward();
        return grads
        
        

class Constant(DObject):
    def __init__(self,owner,value,name=None):
        super().__init__(owner);
        self.owner.n_cte+=1;
        if name is None:
            self.name="cte_"+str(self.owner.n_cte);
        else:
            if name in list(self.owner.cte.keys()):
                self.name="cte_"+str(self.owner.n_cte);
                warnings.warn("Warning the cte with name: "+name+" has been already defined, replacing with "+self.name)
            else:
                self.name=name;
        self.value=value;
        self.gradient=None;
        self.gradient_obj=None;
        # Updating the owner dictionary
        # and assign the object to the dictionary position
        self.owner.cte[self.name]=self;
        self=self.owner.cte[self.name];


    def __repr__(self):
        return "The cte is : name "+self.name+" with value "+str(self.value)

class Variable(DObject):
    def __init__(self,owner,value,name=None):
        super().__init__(owner);
        self.owner.n_var+=1;
        if name is None:
            self.name="var_"+str(self.owner.n_var);
        else:
            if name in list(self.owner.var.keys()):
                self.name="cte_"+str(self.owner.n_var);
                warnings.warn("Warning the variable with name: "+name+" has been already defined, replacing with "+self.name)
            else:
                self.name=name;
        self.value=value;
        # Updating the owner dictionary
        # and assign the object to the dictionary position
        self.gradient=None;
        self.gradient_obj=None;
        self.owner.var[self.name]=self;
        self=self.owner.var[self.name];

    def __repr__(self):
        return "The variable is: name "+self.name+" with value "+str(self.value);


class Operator(DObject):
    def __init__(self,owner,name="ope_"):
        super().__init__(owner);
        self.value=None;
        self.value_obj=None;
        self.inputs=[];
        self.name=name;
        self.gradient=None;
        self.gradient_obj=None;


    def __repr__(self):
        return "The name of the operator is, "+self.name ;

class Add(Operator):
    def __init__(self,owner,a,b,name="add_"):
        super().__init__(owner,name);
        self.owner.n_add+=1;
        self.name="add_"+str(self.owner.n_add)
        self.inputs=[a,b]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name]

    def forward(self,a,b):
        return a+b;

    def forward_obj(self,a,b):
        return a+b;

    def backward(self,a,b,dout):
        return dout,dout;
    
    def backward_obj(self,a,b,dout):
        return dout,dout;

class Sub(Operator):
    def __init__(self,owner,a,b,name="sub_"):
        super().__init__(owner,name);
        self.owner.n_sub+=1;
        self.name="sub_"+str(self.owner.n_sub)
        self.inputs=[a,b]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name]

    def forward(self,a,b):
        return a-b;

    def forward_obj(self,a,b):
        return a-b;

    def backward(self,a,b,dout):
        return dout,-dout;

    def backward_obj(self,a,b,dout):
        return dout,-dout;


class Multiply(Operator):
    def __init__(self,owner,a,b,name="mul_"):
        super().__init__(owner,name);
        self.owner.n_mul+=1;
        self.name="mul_"+str(self.owner.n_mul);
        self.inputs=[a,b]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];


    def forward(self,a,b):
        return a*b;

    def forward_obj(self,a,b):
        return a*b;

    def backward(self,a,b,dout):
        return dout*b,dout*a;

    def backward_obj(self,a,b,dout):
        return dout*b,dout*a

class Power(Operator):
    def __init__(self,owner,a,b,name="pow_"):
        super().__init__(owner,name);
        self.owner.n_pow+=1;
        self.name="pow_"+str(self.owner.n_pow);
        self.inputs=[a,b]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];

    def forward(self,a,b):
        return np.power(a,b)

    def forward_obj(self,a,b):
        return self.owner.set_ope(Power,a,b);

    def backward(self,a,b,dout):
        return dout*b*np.power(a,b-1.0),dout*np.log(a)*np.power(a,b)

    def backward_obj(self,a,b,dout):
        return dout*b*self.owner.set_ope(a,b-1.0),dout*self.owner.set_ope();

class Divide(Operator):
    def __init__(self,owner,a,b,name="div_"):
        super().__init__(owner,name);
        self.owner.n_div+=1;
        self.name="div_"+str(self.owner.n_div);
        self.inputs=[a,b]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];

    def forward(self,a,b):
        return a/b;
    
    def forward_obj(self,a,b):
        return a/b;

    def backward(self,a,b,dout):
        return dout/b,-a*dout/b**2.0;

    def backward_obj(self,a,b,dout):
        return dout/b,-a*dout/b**2.0

class Sine(Operator):
    def __init__(self,owner,a,b=None,name="sine_"):
        super().__init__(owner,name);
        self.owner.n_sine+=1;
        self.name="sine_"+str(self.owner.n_sine);
        self.inputs=[a]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];
    
    def forward(self,a):
        return np.sin(a);

    def backward(self,a,dout):
        return [dout*np.cos(a)];


class Cosine(Operator):
    def __init__(self,owner,a,b=None,name="cos_"):
        super().__init__(owner,name);
        self.owner.n_cos+=1;
        self.name="cos_"+str(self.owner.n_cos);
        self.inputs=[a]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];
    # Here the Cosine
    def forward(self,a):
        return np.cos(a)

    def backward(self,a,dout):
        return [-dout*np.sin(a)]

class LogN(Operator):
    def __init__(self,owner,a,b=None,name="logn_"):
        super().__init__(owner,name);
        self.owner.n_logn+=1;
        self.name="logn_"+str(self.owner.n_logn);
        self.inputs=[a]
        # Add the operator to the owner
        self.owner.ope[self.name]=self;
        self=self.owner.ope[self.name];
    
    def forward(self,a):
        return np.log(a);

    def backward(self,a,dout):
        return [dout/a];


class Tangent(Operator):
    # Here the Tangent
    pass





if __name__=='__main__':
    import sys,os;
    '''
    Just make the autodifferentiation
    with the header that you are looking.
    '''
    #myform = DTrack();
    #x  = myform.set_var(value=5.0,name="x");
    #y  = myform.set_var(value=3.0,name="y") 
    #c1 = myform.set_cte(value=2.0);
    #c2 = myform.set_cte(value=3.0);
    #header = (x**c1+c2*x+c1*(x+c2))/x
    #header = x**2+y**2+(x * y);
    # Maybe the order should be an internal variable
    #myform.set_header(header);
    #order = myform.Topsort();
    #myform.Forward()
    #myform.Backward()
    #myform.PlotGraph();
    # Now get the minimization
    #myform.GetGradObj();
    #print(myform.var["x"].gradient_obj)
    #myform.set_header(myform.var["x"].gradient_obj)
    #myform.Topsort();
    #myform.Forward();

    #grads=myform.Backward(order);
    #print(x.gradient)
    #myform.reset()
    #figure=myform.PlotGraph();
