from numpy import sin, cos, pi
from astropy import coordinates, constants
from astropy import units as u 
from astropy.time import Time
import numpy as np
import lal, lalsimulation

#-----------------------------------------------------COORDIANTE TRANSFORMATION-----------------------------------------------------    
def from_icrs_to_gcrs(icrs_coord):
  if isinstance(icrs_coord,np.ndarray) and gcrs_coord.shape==(3,):
    x,y,z=icrs_coord
    return coordinates.SkyCoord(x,y,z,unit='AU',representation_type='cartesian',frame='icrs').transform_to('gcrs')
  elif icrs_coord.frame is 'icrs':
    return icrs_coord.transform_to('gcrs')
  else :
    raise RuntimeError("1")

def from_gcrs_to_icrs(gcrs_coord):
  if isinstance(gcrs_coord,np.ndarray) and gcrs_coord.shape==(3,):#Add unit section
    x,y,z=gcrs_coord
    return coordinates.SkyCoord(x,y,z,unit='AU',representation_type='cartesian',frame='gcrs').transform_to('icrs')
  elif gcrs_coord.frame is 'gcrs':
    return gcrs_coord.transform_to('icrs')
  else :
    raise RuntimeError("1")

class LISA(object):
  def __init__(self,t_gps,kappa,_lambda_):
    self.t_gps=t_gps
    self.kappa=kappa
    self._lambda_=_lambda_
#-----------------------------------------------------DETECTOR POSITION-----------------------------------------------------
  def get_pos_detector(self,plot=False):
    t=Time(val=self.t_gps,format='gps',scale='utc').to_datetime(timezone=None)
    t_ref = np.array([2034-t.year,t.month/12,t.day/(12*365),t.hour/(12*365*24),t.minute/(12*365*24*60),t.second/(12*365*24*60*60),t.microsecond/(12*365*24*60*60*1e-6)])
    t_ref = np.sum(t_ref,axis=0)
    n=np.array(range(1,4))
    alpha=2.*pi*t_ref/1+self.kappa
    beta_n=(n-1)+2.*pi/3+self._lambda_
    a, L = 1.,.1    #*u.AU
    e = L/(2.*a*np.sqrt(3))

    #       pos[0],pos[1],pos[2] = X, Y, Z for all 3 detectors at one time

    pos = np.array([a*cos(alpha) + a*e*(sin(alpha)*cos(alpha)*sin(beta_n)-(1 + sin(alpha)**2)*cos(beta_n)),
                    a*sin(alpha) + a*e*(sin(alpha)*cos(alpha)*sin(beta_n)-(1 + cos(alpha)**2)*sin(beta_n)),
                    -np.sqrt(3)*a*e*cos(alpha - beta_n)])
    
    if plot:

      from mpl_toolkits.mplot3d import Axes3D
      import matplotlib.pyplot as plt 
      
      #ax.scatter(pos[0],pos[1],pos[2], marker='o')# X,Y,Z at current time 
      t=np.arange(0,10,.1)
      x_1=a*cos(t) + a*e*(sin(t)*cos(t)*sin(beta_n[0])-(1 + sin(t)**2)*cos(beta_n[0]))
      y_1=a*sin(t) + a*e*(sin(t)*cos(t)*sin(beta_n[0])-(1 + cos(t)**2)*sin(beta_n[0]))
      z_1=-np.sqrt(3)*a*e*cos(t - beta_n[0])
      x_2=a*cos(t) + a*e*(sin(t)*cos(t)*sin(beta_n[1])-(1 + sin(t)**2)*cos(beta_n[1]))
      y_2=a*sin(t) + a*e*(sin(t)*cos(t)*sin(beta_n[1])-(1 + cos(t)**2)*sin(beta_n[1]))
      z_2=-np.sqrt(3)*a*e*cos(t - beta_n[1])
      x_3=a*cos(t) + a*e*(sin(t)*cos(t)*sin(beta_n[2])-(1 + sin(t)**2)*cos(beta_n[2]))
      y_3=a*sin(t) + a*e*(sin(t)*cos(t)*sin(beta_n[2])-(1 + cos(t)**2)*sin(beta_n[2]))
      z_3=-np.sqrt(3)*a*e*cos(t - beta_n[2])
      rand_pt=from_gcrs_to_icrs(np.zeros(3))
      t=Time(val=self.t_gps,format='gps',scale='utc')
      sun=coordinates.get_sun(t).transform_to(frame='icrs')
      sun.representation_type, rand_pt.representation_type='cartesian', 'cartesian'
      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      ax.scatter(np.float32(rand_pt.x),np.float32(rand_pt.y),np.float32(rand_pt.z),marker=',')
      ax.scatter(np.float32(sun.x),np.float32(sun.y),np.float32(sun.z), marker='h')
      ax.scatter(x_1,y_1,z_1, marker='o')
      ax.scatter(x_2,y_2,z_2, marker='+')
      ax.scatter(x_3,y_3,z_3, marker='*')

    coord_ICRS=coordinates.SkyCoord(pos[0],pos[1],pos[2],unit=u.AU,representation_type='cartesian',frame='icrs')
    return coord_ICRS
    #return np.array([x,y,z])
    
#--------------------------------------------------DISTANCE FROM DETECTOR-------------------------------------------

  def light_travel_time_to_detector(self,det,ref_time):
    if isinstance(det,str):              #if ref_time is None:
      det_loc=from_gcrs_to_icrs(lalsimulation.DetectorPrefixToLALDetector('H1').location*6.6846e-12)
      det_loc.representation_type='cartesian'
      _a_=np.array([np.float32(det_loc.x),np.float32(det_loc.y),np.float32(det_loc.z)])
      L_pos=LISA(ref_time,self.kappa,self._lambda_).get_pos_detector()
      _b_=np.array([np.float32(L_pos.x),np.float32(L_pos.y),np.float32(L_pos.z)])
      d=_a_-_b_
      return d.dot(d*0.5)/constants.c.value
      
#--------------------------------------------------DISTANCE FROM LOCATION-------------------------------------------

  def light_time_delay_from_location(self,ref_time,other_location): #similar to time_delay_from_location from pycbc.detector
    L_pos=LISA(ref_time,self.kappa,self._lambda_).get_pos_detector()             # UNIT AND Coordinate system and time
    _b_=np.array([np.float32(L_pos.x),np.float32(L_pos.y),np.float32(L_pos.z)])
    if isinstance(other_location,np.ndarray):
      _a_=other_location
    elif isinstance(other_location,coordinates.SkyCoord):
      _a_=np.array([np.float32(other_location.x),np.float32(other_location.y),np.float32(other_location.z)])
    d=_a_-_b_
    return d.dot(d*0.5)/constants.c.value
