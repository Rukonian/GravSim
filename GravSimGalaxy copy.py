import numpy as np
import cython
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')

R = 8.3145

k = 1.380649E-23

pi = np.pi

G = 6.67430E-11

class OctNode:
    def __init__(self, center, size, masses, points, ids, leaves=[]):
        self.center=center #center of the node
        self.size=size #maximum size length of the node
        self.children=[] #empty list, nodes do not start with child nodes
    
        Npoints=len(points)

        if Npoints==1: #if there is only one point, we create a node
            leaves.append(self)
            print(leaves)
            print(' ')
            self.COM=points[0]
            self.mass=masses[0]
            self.id=ids[0]
            self.g=np.zeros(3)
        else: #if a node has 2+ points, create child nodes
            self.GenerateChildren(points, masses, ids, leaves)
            com_total=np.zeros(3) #running com moments
            m_total=0 #running total node mass
            for c in self.children:
                m, com=c.mass, c.COM
                m_total+=m
                com_total+=com*m #sum moments in node
            self.mass=m_total
            self.COM=com_total/self.mass

    def GenerateChildren(self, points, masses, ids, leaves):
        octant_index=(points>self.center) #determines point octant
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    in_octant=np.all(octant_index==np.bool_([i,j,k]), axis=1)
                    if not np.any(in_octant):
                        continue #dont make a node without particles
                    dx=0.5*self.size*(np.asarray([i,j,k])-0.5)
                    self.children.append(OctNode(self.center+dx, self.size/2, masses[in_octant],
                                                 points[in_octant], ids[in_octant], leaves))

def TreeWalk(node, node0, thetamax=0.7, G=G):
    dx=node.COM-node0.COM #vector between node COMs
    r=np.sqrt(np.sum(dx**2)) #distance for thetamax calc
    if r>0:
        if (len(node.children)==0) or (node.size/r<thetamax): #if the node has one particle or theta small enough
            node0.g+=G*node.mass*dx/(r**3)
        else: #otherwise, recursively go through child nodes until all nodes done
            for c in node.children: 
                TreeWalk(c, node0, thetamax, G)

def GravAccel(points, masses, thetamax=0.7, G=G):
    center=(np.max(points, axis=0)+np.min(points, axis=0))/2
    topsize=np.max(np.max(points, axis=0)-np.min(points, axis=0))
    leaves=[]
    topnode=OctNode(center, topsize, masses, points, np.arange(len(points)), leaves)

    accel=np.empty_like(points)
    for i, leaf in enumerate(leaves):
        TreeWalk(topnode, leaf, thetamax, G)
        accel[leaf.id]=leaf.g
    return accel


class particle: #Defines the particle class, meant for simple gravitational and EM simulations
    def __init__(self, t, x, y, z, vx, vy, vz, mass=2E30, radius=1.73E-10, charge=0):
        self.mass=mass
        self.radius=radius
        self.charge=charge
        self.t=t
        self.x=x
        self.y=y
        self.z=z
        self.vx=vx
        self.vy=vy
        self.vz=vz
        
    def __getitem__(self, item):
        if item=='x':
            return (self.t, self.x, self.y, self.z)
        if item=='u':
            return (self.t, self.vx, self.vy, self.vz)
        
class System:
    def __init__(self, n, T, size, dt=0.5): #number of particles, kinetic temperature in kelvin, bounding cub radius in meters, timestep in seconds
        self.n=n
        self.T=T
        self.size=size
        self.dt=dt
        self.bodies=[] #stores all particle objects in sim
        self.points=[] #stores particles points as list of np.arrays
        self.masses=[] #stores particles masses as list of np.arrays
        self.velocities=[]
        
    def add_body(self): #Adds particles to the simulation
        for i in range(self.n):
            self.bodies.append(particle(0, np.random.uniform(-self.size, self.size), np.random.uniform(-self.size, self.size), np.random.uniform(-self.size, self.size),  #x, y, z coords
                                        (1 if np.random.random() >=0.5 else -1)*np.random.normal(np.sqrt((8*R*self.T)/(pi*29)), np.sqrt(self.T*(k)/(4.82E-26))), (1 if np.random.random() >=0.5 else -1)*np.random.normal(np.sqrt((8*R*self.T)/(pi*29)), np.sqrt(self.T*(k)/(4.82E-26))),  #vx, vy, vz distribution
                                                                                                                                                  (1 if np.random.random() >=0.5 else -1)*np.random.normal(np.sqrt((8*R*self.T)/(pi*29)), np.sqrt(self.T*(k)/(4.82E-26)))))
            
    def advance(self, duration, thetamax=0.5):
        time=0
        for body in self.bodies:
            self.points.append([body['x'][1], body['x'][2], body['x'][3]])
            self.velocities.append([body['u'][1], body['u'][2], body['u'][3]])
            self.masses.append(body.mass)
            body.t+=self.dt
        points=np.asarray(self.points)
        masses=np.asarray(self.masses)
        velocities=np.asarray(self.velocities)
        #df0=pd.DataFrame({'(x,y,z)': points, 'mass': masses, 'velocity': velocities})
        #pd.to_csv(r'C:\Users\Andym\OneDrive\Documents\PythonGravSimData\t0')
        while time<duration:
            a=GravAccel(points, masses, thetamax)

            #print(np.sqrt(np.sum(points**2, axis=1)))
            
            velocities+=(a)*self.dt #Calculating new velocity vector based on force field and current velocity vector
            points+=velocities*self.dt #moving particle according to updated velocity vector

            time+=1

            #x=points[:,0]
            #y=points[:,1]
            #z=points[:,2]
            #
            #fig = plt.figure()
            #plt.figure(figsize=(12,12))
            #ax = fig.add_subplot(projection='3d')
            #ax.scatter(x, y, z, c='blue', alpha=0.8, s=1, marker='.')
            #
            ## Set an equal aspect ratio
            #ax.set_aspect('equal')
            #fig.savefig(r'C:\Users\Andym\OneDrive\Documents\PythonGravSimData\Python\Galaxy{}.png'.format(time), dpi=200, format='png', bbox_inches='tight')
            #plt.close('all')
            #if time%(5)==0:
                #print(time/(duration/self.dt))
            
            #df=pd.DataFrame({'(x,y,z)': points, 'mass': masses, 'velocity': velocities})
            #df.to_csv(r'C:\Users\Andym\OneDrive\Documents\PythonGravSimData\t{}'.format(str(time).replace('.','')))

sim=System(5, 5, 5, dt=1)

sim.add_body()

sim.advance(5)
