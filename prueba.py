class graphic:
    def __init__(self):
        self.dic_perros={};
    def __repr__(self):
        return "grafico"


def sin(perro):
    print(perro.nombre)

class perro:
    def __init__(self,owner,nombre):
        self.owner=owner;
        self.nombre=nombre;
        print(self.owner)
        self.owner.dic_perros[self.nombre]=self;
        self=self.owner.dic_perros[self.nombre];

    def __repr__(self):
        return "Name: "+self.nombre


migrafico=graphic()
a=perro(migrafico,"perro")
a1=perro(migrafico,"perro1")
print(migrafico.dic_perros)
a1.nombre="perro2"
a2=a1
a2.nombre="perro3"
print(migrafico.dic_perros)
sin(a)
print("{:.2f}".format(2.13333))


def function():
    x=1
    def _printea():
        print(x)
    _printea()


function()
