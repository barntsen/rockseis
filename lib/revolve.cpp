#include "revolve.h"

namespace rockseis {
// constructor
//
template<typename T>
Revolve<T>::Revolve(unsigned int _nt, unsigned int _snaps)
{
    steps = _nt;
    snaps = _snaps;
    if (snaps > steps || snaps < 2) {
        rs_warning("snaps invalid - setting to default optimal value.");
        snaps = adjust(steps);
    }
    
    chp = (T *) calloc(1,1);
    Fc = std::make_shared<File>(); 

	/* Neccessary for first call to revolve */
	capo=0;
	fine = steps;
	check = -1;    
}

template <typename T>
int Revolve<T>::numforw(int steps, int snaps)
/*<Number of forward computing steps needed>*/
{
    int reps, range, num;

    if (snaps < 1) 
    {
        fprintf(stderr," error occurs in numforw: snaps < 1\n"); 
        return -1;
    }
    if (snaps > CHECKUP)
    {
        fprintf(stderr," number of snaps=%d exceeds CHECKUP \n",snaps);
        fprintf(stderr," redefine 'CHECKUP' \n");
        return -1;
    }
    reps = 0;
    range = 1;
    while(range < steps)
    { 
        reps += 1;
        range = range*(reps + snaps)/reps; 
    }
    if (reps > REPSUP)
    {
        fprintf(stderr," number of reps=%d exceeds REPSUP \n",reps);
        fprintf(stderr," redefine 'REPSUP' \n");
        return -1;
    }
    num = reps * steps - range*reps/(snaps+1);
    return num;
}

template <typename T>
double Revolve<T>::expense(int steps, int snaps)
 {
  double ratio;

  if (snaps < 1)
   {
    fprintf(stderr," error occurs in expense: snaps < 0\n");
    return -1;
   }
  if (steps < 1)
   {
    fprintf(stderr," error occurs in expense: steps < 0\n");
    return -1;
   }
  ratio = ((double) numforw(steps,snaps));
  if (ratio == -1)
    return -1;
  ratio = ratio/steps;
  return ratio;
 }

template <typename T>
int Revolve<T>::maxrange(int ss, int tt)
/*<returns the integer (SNAPS+REPS)!/(SNAPS!REPS!) provided SNAPS >=0, REPS >= 0>*/
 {
  int i, ires;
  double res = 1.0;

  if((tt<0) || (ss<0))
   {
    fprintf(stderr,"error in MAXRANGE: negative parameter");
    return -1;
   }
  for(i=1; i<= tt; i++)
   {  
    res *= (ss + i);
    res /= i;
    if (res > MAXINT)
     { 
      ires=MAXINT;
      fprintf(stderr,"warning from MAXRANGE: returned maximal integer %d\n",ires);
      return ires;
     }
   }
  ires = res;
  return ires;
 }

template <typename T>
int Revolve<T>::adjust(int steps)
/*<It can be used to determine a value of SNAPS so that the increase in spatial complexity equals approximately the increase in temporal complexity.>*/
 {
  int snaps, s, reps;

  snaps = 1;
  reps = 1;
  s = 0;
  while( maxrange(snaps+s, reps+s) > steps ) 
    s--;
  while( maxrange(snaps+s, reps+s) < steps ) 
    s++;
  snaps += s;
  reps += s ;
  s = -1;
  while( maxrange(snaps,reps) >= steps )
   {
    if (snaps > reps)
     { 
      snaps -= 1; 
      s = 0; 
     }
    else 
     { 
      reps -= 1; 
      s = 1; 
     }
   }
  if ( s == 0 ) 
    snaps += 1 ;
  if ( s == 1 ) 
    reps += 1;
  return snaps;
 }

template <typename T>
revolve_action Revolve<T>::revolve()
/*<Control for the forward and reverse state flows>*/
 {
  int ds, oldcapo, num, bino1, bino2, bino3, bino4, bino5;
    /* (capo,fine) is the time range currently under consideration */
    /* ch[j] is the number of the state that is stored in checkpoint j */

  numbers.commands += 1;
  if ((check < -1) || (capo > fine)) 
    return error; 
  if ((check == -1) && (capo < fine))
   {
    if (check == -1) 
      turn = 0;   /* initialization of turn counter */
    *ch = capo-1;   
   }
  switch(fine-capo)
   { 
    case 0:   /* reduce capo to previous checkpoint, unless done  */
      if(check == -1 || capo==*ch )
       { 
        check -= 1;
        if (info > 0)  
         { 
          fprintf(stderr," \n advances: %5d",numbers.advances);
          fprintf(stderr," \n takeshots: %4d",numbers.takeshots);
          fprintf(stderr," \n commands: %5d \n",numbers.commands);
         } 
        return terminate;
       }
      else
       { 
        capo = ch[check];
        oldfine = fine;
        return restore;
       } 
    case 1:  /* (possibly first) combined forward/reverse step */ 
      fine -= 1;
      if(check >= 0 && ch[check] == capo) 
        check -= 1; 
      if(turn == 0)
       {
        turn = 1;
        oldfine = fine;
        return firsturn; 
       }
      else 
       { 
        oldfine = fine;
        return youturn; 
       }
    default:         
      if(check == -1 || ch[check] != capo) 
       { 
        check += 1 ; 
        if(check >= CHECKUP)
	 { 
	  info = 10;
          return error;
	 }
        if(check+1 > snaps)
         {
          info = 11;
          return error;
         }
        ch[check] = capo;
        if (check == 0) 
         {
          numbers.advances = 0;
          numbers.takeshots = 0;
          numbers.commands = 1;
          oldsnaps = snaps;
          if (snaps > CHECKUP)
           {
            info = 14;
            return error;
           }
          if (info > 0) 
           {
            num = numforw(fine-capo,snaps);
            if (num == -1) 
             {
              info = 12;
              return error;
             }
            fprintf(stderr," prediction of needed forward steps: %8d => \n",num);
            fprintf(stderr," slowdown factor: %8.4f \n\n",((double) num)/(fine-capo));
           }
         }
        numbers.takeshots += 1;
        oldfine = fine;
        return takeshot; 
       }
      else
       { 
        if ((oldfine < fine) && (snaps == check+1))
         { 
          info = 13;
          return error;
         } 
        oldcapo = capo;
        ds = snaps - check;
        if (ds < 1)
         {
          info = 11;
          return error;
         }
        reps = 0;
        range = 1;
        while(range < fine - capo) 
         { 
          reps += 1;
          range = range*(reps + ds)/reps; 
         }
        if (reps > REPSUP) 
         { 
          info = 15;
          return error;
         }
        if (snaps != oldsnaps)
         { 
           if (snaps > CHECKUP)
            {
             info = 14;
             return error;
            }
         }
        bino1 = range*reps/(ds+reps);
        bino2 = (ds > 1) ? bino1*ds/(ds+reps-1) : 1;
        if (ds == 1)
          bino3 = 0;
        else
          bino3 = (ds > 2) ? bino2*(ds-1)/(ds+reps-2) : 1;
        bino4 = bino2*(reps-1)/ds;
        if (ds < 3)
          bino5 = 0;
        else
          bino5 = (ds > 3) ? bino3*(ds-2)/reps : 1;
     
        if (fine-capo <= bino1 + bino3)
          capo = capo+bino4;
        else 
         {
          if (fine-capo >= range - bino5) 
            capo = capo + bino1; 
          else 
             capo = fine-bino2-bino3;
         }
        if (capo == oldcapo) 
          capo = oldcapo+1;  
        numbers.advances = numbers.advances + capo - oldcapo; 
        oldfine = fine;
        return advance;
       }          
    }
 } 

// destructor
template<typename T>
Revolve<T>::~Revolve(){
    free(chp);
}

// =============== INITIALIZING TEMPLATE STRUCTS =============== //
template class Revolve<float>;
template class Revolve<double>;

}


