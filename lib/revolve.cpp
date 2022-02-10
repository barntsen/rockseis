#include "revolve.h"

// Regex s/memcpy(\(.\{-}\) \(.\{-},\)/memcpy(\2 \1/
namespace rockseis {
// constructor
//
template<typename T>
Revolve<T>::Revolve(unsigned int _nt, unsigned int _snaps, bool _incore)
{
   steps = _nt;
   snaps = _snaps;
   if (snaps > steps || snaps < 2) {
      rs_warning("snaps invalid - setting to default optimal value.");
      snaps = adjust(steps);
   }

   Fc = std::make_shared<File>(); 

   /* Neccessary for first call to revolve */
   capo=0;
   fine = steps;
   check = -1;    

   info = 0;
   incore = _incore;
   open = false;
   allocated = false;
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

template<typename T>
void Revolve<T>::createCheck(std::string _filename, char flag)
{
   switch(flag){
      case 'w':
         if(!incore){
            this->filename = _filename;
            if(!this->filename.empty()){
               this->Fc->output(this->filename);
               this->open = true;
               this->Fc->setN(1, this->checksize);
               this->Fc->setN(2, this->snaps);
               this->Fc->setData_format(sizeof(T));
               this->Fc->setType(rockseis::CHECKPOINT);
               this->Fc->writeHeader();
               this->Fc->seekp(this->Fc->getStartofdata());
            }else{
               rs_error("Revolve::createCheck: Filename not set.");
            }
         }else{
            this->checkpoints = (T *) calloc(this->checksize*this->snaps, sizeof(T));
            if(checkpoints == NULL) rs_error("Revolve::createCheck: Failed to allocate memory for checkpoints.");
            this->allocated = true;
         }
         break;
      case 'r':
         if(this->open) rs_error("Checkpoint cannot be opened two times.");
         this->filename = _filename;
         if(this->Fc->input(this->filename) == FILE_ERR)
         {
            rs_error("Revolve::createCheck: Error opening checkpoint file for reading.");
         }
         if(this->checksize != this->Fc->getN(1)) rs_error("Revolve::createCheck: Mismatch in size of checkpoints");
         if(this->snaps != this->Fc->getN(2)) rs_error("Revolve::createCheck: Mismatch in number of checkpoints");
         if(sizeof(T) != this->Fc->getData_format()) rs_error("Revolve::createCheck: Mismatch in precision of checkpoints");
         if(this->Fc->getType() != rockseis::CHECKPOINT) rs_error("Revolve::createCheck: Mismatch in file type");
         if(this->incore){
            if(this->allocated) free(this->checkpoints);
            this->checkpoints = (T *) calloc(this->checksize*this->snaps, sizeof(T));
            if(checkpoints == NULL) rs_error("Revolve::createCheck: Failed to allocate memory for checkpoints.");
         }else{
            this->open = true;
         }
         break;
      case 'a':
         if(!this->incore){
            if(this->open) rs_error("Checkpoint cannot be opened two times.");
            this->filename = _filename;
            if(this->Fc->append(this->filename) == FILE_ERR)
            {
               rs_error("Revolve::createCheck: Error opening checkpoint file for reading and writting.");
            }
            if(this->checksize != this->Fc->getN(1)) rs_error("Revolve::createCheck: Mismatch in size of checkpoints");
            if(this->snaps != this->Fc->getN(2)) rs_error("Revolve::createCheck: Mismatch in number of checkpoints");
            if(sizeof(T) != this->Fc->getData_format()) rs_error("Revolve::createCheck: Mismatch in precision of checkpoints");
            if(this->Fc->getType() != rockseis::CHECKPOINT) rs_error("Revolve::createCheck: Mismatch in file type");
            this->open = true;
         }
         break;
      default: 
         rs_error("Revolve::createCheck: Invalid flag.");
   }
}


template<typename T>
void Revolve<T>::openCheck(std::string _filename, std::shared_ptr<WavesAcoustic2D<T>> waves, char flag)
{
   long long nx_pml, nz_pml;
   long long lpml;

   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   this->checksize = 3*nz_pml*nx_pml;
   if((waves->getPml())->getApplypml(0)){
      this->checksize += 2*nz_pml*lpml;
   }
   if((waves->getPml())->getApplypml(1)){
      this->checksize += 2*nz_pml*lpml;
   }
   if((waves->getPml())->getApplypml(4)){
      this->checksize += 2*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(5)){
      this->checksize += 2*nx_pml*lpml;
   }

   this->createCheck(_filename, flag); 
}

template<typename T>
void Revolve<T>::openCheck(std::string _filename, std::shared_ptr<WavesAcoustic3D<T>> waves, char flag)
{
   long long nx_pml, ny_pml, nz_pml;
   long long lpml;

   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   this->checksize = 4*nz_pml*nx_pml*ny_pml;
   if((waves->getPml())->getApplypml(0)){
      this->checksize += 2*nz_pml*ny_pml*lpml;
   }
   if((waves->getPml())->getApplypml(1)){
      this->checksize += 2*nz_pml*ny_pml*lpml;
   }
   if((waves->getPml())->getApplypml(2)){
      this->checksize += 2*nz_pml*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(3)){
      this->checksize += 2*nz_pml*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(4)){
      this->checksize += 2*nx_pml*ny_pml*lpml;
   }
   if((waves->getPml())->getApplypml(5)){
      this->checksize += 2*nx_pml*ny_pml*lpml;
   }

   if(sizeof(T)*this->checksize > LIMIT_WARNING)
   {
      rs_warning("Revolve<T>::openCheck::Required memory per core exceeds 2 Gb: ", std::to_string(4.*this->checksize/1024./1024./1024.));
   }
   this->createCheck(_filename, flag); 
}

template<typename T>
void Revolve<T>::openCheck(std::string _filename, std::shared_ptr<WavesElastic2D<T>> waves, char flag)
{
   long long nx_pml, nz_pml;
   long long lpml;

   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   this->checksize = 5*nz_pml*nx_pml;

   if((waves->getPml())->getApplypml(0)){
      this->checksize += 4*nz_pml*lpml;
   }
   if((waves->getPml())->getApplypml(1)){
      this->checksize += 4*nz_pml*lpml;
   }
   if((waves->getPml())->getApplypml(4)){
      this->checksize += 4*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(5)){
      this->checksize += 4*nx_pml*lpml;
   }
   this->createCheck(_filename, flag); 
}

template<typename T>
void Revolve<T>::openCheck(std::string _filename, std::shared_ptr<WavesElastic2D_DS<T>> waves, char flag)
{
   long long nx_pml, nz_pml;
   long long lpml;

   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   this->checksize = 4*nz_pml*nx_pml + 8*nz_pml*lpml + 8*nx_pml*lpml;
   this->createCheck(_filename, flag); 
}

template<typename T>
void Revolve<T>::openCheck(std::string _filename, std::shared_ptr<WavesElastic3D<T>> waves, char flag)
{
   long long nx_pml, ny_pml, nz_pml;
   long long lpml;

   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   this->checksize = 9*nz_pml*nx_pml*ny_pml;
   if((waves->getPml())->getApplypml(0)){
      this->checksize += 6*nz_pml*ny_pml*lpml;
   }
   if((waves->getPml())->getApplypml(1)){
      this->checksize += 6*nz_pml*ny_pml*lpml;
   }
   if((waves->getPml())->getApplypml(2)){
      this->checksize += 6*nz_pml*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(3)){
      this->checksize += 6*nz_pml*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(4)){
      this->checksize += 6*nx_pml*ny_pml*lpml;
   }
   if((waves->getPml())->getApplypml(5)){
      this->checksize += 6*nx_pml*ny_pml*lpml;
   }

   if(sizeof(T)*this->checksize > LIMIT_WARNING)
   {
      rs_warning("Revolve<T>::openCheck::Required memory per core exceeds 2 Gb: ", std::to_string(4.*this->checksize/1024./1024./1024.));
   }

   this->createCheck(_filename, flag); 
}

template<typename T>
void Revolve<T>::openCheck(std::string _filename, std::shared_ptr<WavesElastic3D_DS<T>> waves, char flag)
{
   long long nx_pml, ny_pml, nz_pml;
   long long lpml;

   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   this->checksize = 6*nz_pml*nx_pml*ny_pml + 12*nz_pml*ny_pml*lpml + 12*nz_pml*nx_pml*lpml + 12*nx_pml*ny_pml*lpml;
   if(sizeof(T)*this->checksize > LIMIT_WARNING)
   {
      rs_warning("Revolve<T>::openCheck::Required memory per core exceeds 2 Gb: ", std::to_string(4.*this->checksize/1024./1024./1024.));
   }
   this->createCheck(_filename, flag); 
}

template<typename T>
void Revolve<T>::openCheck(std::string _filename, std::shared_ptr<WavesViscoelastic2D<T>> waves, char flag)
{
   long long nx_pml, nz_pml;
   long long lpml;

   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   this->checksize = 8*nz_pml*nx_pml;
   if((waves->getPml())->getApplypml(0)){
      this->checksize += 4*nz_pml*lpml;
   }
   if((waves->getPml())->getApplypml(1)){
      this->checksize += 4*nz_pml*lpml;
   }
   if((waves->getPml())->getApplypml(4)){
      this->checksize += 4*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(5)){
      this->checksize += 4*nx_pml*lpml;
   }

   this->createCheck(_filename, flag); 
}

template<typename T>
void Revolve<T>::openCheck(std::string _filename, std::shared_ptr<WavesVti2D<T>> waves, char flag)
{
   long long nx_pml, nz_pml;
   long long lpml;

   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   this->checksize = 5*nz_pml*nx_pml;

   if((waves->getPml())->getApplypml(0)){
      this->checksize += 4*nz_pml*lpml;
   }
   if((waves->getPml())->getApplypml(1)){
      this->checksize += 4*nz_pml*lpml;
   }
   if((waves->getPml())->getApplypml(4)){
      this->checksize += 4*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(5)){
      this->checksize += 4*nx_pml*lpml;
   }
   this->createCheck(_filename, flag); 
}

template<typename T>
void Revolve<T>::openCheck(std::string _filename, std::shared_ptr<WavesOrtho3D<T>> waves, char flag)
{
   long long nx_pml, ny_pml, nz_pml;
   long long lpml;

   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   this->checksize = 9*nz_pml*nx_pml*ny_pml;
   if((waves->getPml())->getApplypml(0)){
      this->checksize += 6*nz_pml*ny_pml*lpml;
   }
   if((waves->getPml())->getApplypml(1)){
      this->checksize += 6*nz_pml*ny_pml*lpml;
   }
   if((waves->getPml())->getApplypml(2)){
      this->checksize += 6*nz_pml*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(3)){
      this->checksize += 6*nz_pml*nx_pml*lpml;
   }
   if((waves->getPml())->getApplypml(4)){
      this->checksize += 6*nx_pml*ny_pml*lpml;
   }
   if((waves->getPml())->getApplypml(5)){
      this->checksize += 6*nx_pml*ny_pml*lpml;
   }

   if(sizeof(T)*this->checksize > LIMIT_WARNING)
   {
      rs_warning("Revolve<T>::openCheck::Required memory per core exceeds 2 Gb: ", std::to_string(4.*this->checksize/1024./1024./1024.));
   }

   this->createCheck(_filename, flag); 
}


template<typename T>
void Revolve<T>::closeCheck()
{
   if(this->open){
      Fc->close();
      this->open = false;
   }
}

template<typename T>
void Revolve<T>::removeCheck()
{
   if(!incore){
      if(this->open) {
         this->Fc->close();
         this->open = false;
      }
      if(!this->filename.empty()){
         if( remove( filename.c_str() ) != 0 ){
            rs_error( "Snapshot::removeSnap: Error deleting file: ", filename);
         }
      }
   }
}


template<typename T>
void Revolve<T>::readCheck(std::shared_ptr<WavesAcoustic2D<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *P = waves->getP();
   T *Vx = waves->getVx();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlAcoustic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::readCheck: checkpoint array is not allocated.");;
      memcpy(P, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Vx, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Vz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;



      if(Pml->getApplypml(0)){
         memcpy(Pml->P_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vxx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(1)){
         memcpy(Pml->P_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vxx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(4)){
         memcpy(Pml->P_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vzz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;
      }

      if(Pml->getApplypml(5)){
         memcpy(Pml->P_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vzz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->read(P, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->read(Vx, nz_pml*nx_pml);
         this->Fc->read(Vz, nz_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->read(Pml->P_left, nz_pml*lpml);
            this->Fc->read(Pml->Vxx_left, nz_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->read(Pml->P_right, nz_pml*lpml);
            this->Fc->read(Pml->Vxx_right, nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->read(Pml->P_top, nx_pml*lpml);
            this->Fc->read(Pml->Vzz_top, nx_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->read(Pml->P_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Vzz_bottom, nx_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::readCheck: Error reading checkpoints from file.");
      }else{
         rs_error("Revolve::readCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::writeCheck(std::shared_ptr<WavesAcoustic2D<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *P = waves->getP();
   T *Vx = waves->getVx();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlAcoustic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::writeCheck: checkpoint array is not allocated.");;
      memcpy(this->checkpoints+pos, P, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Vx, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Vz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      if(Pml->getApplypml(0)){
         memcpy(this->checkpoints+pos, Pml->P_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(1)){
         memcpy(this->checkpoints+pos, Pml->P_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(4)){
         memcpy(this->checkpoints+pos, Pml->P_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;
      }

      if(Pml->getApplypml(5)){
         memcpy(this->checkpoints+pos, Pml->P_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzz_bottom, nx_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->write(P, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->write(Vx, nz_pml*nx_pml);
         this->Fc->write(Vz, nz_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->write(Pml->P_left, nz_pml*lpml);
            this->Fc->write(Pml->Vxx_left, nz_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->write(Pml->P_right, nz_pml*lpml);
            this->Fc->write(Pml->Vxx_right, nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->write(Pml->P_top, nx_pml*lpml);
            this->Fc->write(Pml->Vzz_top, nx_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->write(Pml->P_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Vzz_bottom, nx_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::writeCheck: Error writing checkpoints to file.");
      }else{
         rs_error("Revolve::writeCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::readCheck(std::shared_ptr<WavesAcoustic3D<T>> waves)
{

   size_t nx_pml, ny_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *P = waves->getP();
   T *Vx = waves->getVx();
   T *Vy = waves->getVy();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlAcoustic3D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::readCheck: checkpoint array is not allocated.");
      memcpy(P, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Vx, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Vy, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Vz, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      if(Pml->getApplypml(0)){
         memcpy(Pml->P_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vxx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(1)){
         memcpy(Pml->P_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vxx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(2)){
         memcpy(Pml->P_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Vyy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(3)){
         memcpy(Pml->P_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Vyy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
      }
      if(Pml->getApplypml(4)){
         memcpy(Pml->P_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(5)){
         memcpy(Pml->P_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->read(P,nz_pml*ny_pml*nx_pml, pos*sizeof(T));
         this->Fc->read(Vx,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Vy,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Vz,nz_pml*ny_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->read(Pml->P_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxx_left,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->read(Pml->P_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxx_right,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(2)){
            this->Fc->read(Pml->P_front,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Vyy_front,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(3)){
            this->Fc->read(Pml->P_back,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Vyy_back,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->read(Pml->P_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzz_top,nx_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->read(Pml->P_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzz_bottom,nx_pml*ny_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::readCheck: Error reading checkpoints from file.");
      }else{
         rs_error("Revolve::readCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::writeCheck(std::shared_ptr<WavesAcoustic3D<T>> waves)
{

   size_t nx_pml, ny_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *P = waves->getP();
   T *Vx = waves->getVx();
   T *Vy = waves->getVy();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlAcoustic3D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::writeCheck: checkpoint array is not allocated.");
      memcpy(this->checkpoints+pos, P, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Vx, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Vy, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Vz, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      if(Pml->getApplypml(0)){
         memcpy(this->checkpoints+pos, Pml->P_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(1)){
         memcpy(this->checkpoints+pos, Pml->P_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(2)){
         memcpy(this->checkpoints+pos, Pml->P_front, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyy_front, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(3)){
         memcpy(this->checkpoints+pos, Pml->P_back, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyy_back, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
      }
      if(Pml->getApplypml(4)){
         memcpy(this->checkpoints+pos, Pml->P_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(5)){
         memcpy(this->checkpoints+pos, Pml->P_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->write(P,nz_pml*ny_pml*nx_pml, pos*sizeof(T));
         this->Fc->write(Vx,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Vy,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Vz,nz_pml*ny_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->write(Pml->P_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxx_left,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->write(Pml->P_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxx_right,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(2)){
            this->Fc->write(Pml->P_front,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Vyy_front,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(3)){
            this->Fc->write(Pml->P_back,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Vyy_back,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->write(Pml->P_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzz_top,nx_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->write(Pml->P_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzz_bottom,nx_pml*ny_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::writeCheck: Error writting checkpoints to file.");
      }else{
         rs_error("Revolve::writeCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::readCheck(std::shared_ptr<WavesElastic2D<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Szz = waves->getSzz();
   T *Sxz = waves->getSxz();
   T *Vx = waves->getVx();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::readCheck: checkpoint array is not allocated.");;
      memcpy(Sxx, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Szz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Sxz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Vx, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Vz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      if(Pml->getApplypml(0)){
         memcpy(Pml->Sxx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Sxzx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vxx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vzx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(1)){
         memcpy(Pml->Sxx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Sxzx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vxx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vzx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(4)){
         memcpy(Pml->Szz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Sxzz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vzz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vxz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;
      }

      if(Pml->getApplypml(5)){
         memcpy(Pml->Szz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Sxzz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vzz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vxz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->read(Sxx, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->read(Szz, nz_pml*nx_pml);
         this->Fc->read(Sxz, nz_pml*nx_pml);
         this->Fc->read(Vx, nz_pml*nx_pml);
         this->Fc->read(Vz, nz_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->read(Pml->Sxx_left, nz_pml*lpml);
            this->Fc->read(Pml->Sxzx_left, nz_pml*lpml);
            this->Fc->read(Pml->Vxx_left, nz_pml*lpml);
            this->Fc->read(Pml->Vzx_left, nz_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->read(Pml->Sxx_right, nz_pml*lpml);
            this->Fc->read(Pml->Sxzx_right, nz_pml*lpml);
            this->Fc->read(Pml->Vxx_right, nz_pml*lpml);
            this->Fc->read(Pml->Vzx_right, nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->read(Pml->Szz_top, nx_pml*lpml);
            this->Fc->read(Pml->Sxzz_top, nx_pml*lpml);
            this->Fc->read(Pml->Vzz_top, nx_pml*lpml);
            this->Fc->read(Pml->Vxz_top, nx_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->read(Pml->Szz_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Sxzz_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Vzz_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Vxz_bottom, nx_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::readCheck: Error reading checkpoints from file.");
      }else{
         rs_error("Revolve::readCheck: File is closed.");

      }
   }
}

   template<typename T>
void Revolve<T>::writeCheck(std::shared_ptr<WavesElastic2D<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Szz = waves->getSzz();
   T *Sxz = waves->getSxz();
   T *Vx = waves->getVx();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::writeCheck: checkpoint array is not allocated.");;
      memcpy(this->checkpoints+pos, Sxx, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Szz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Sxz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Vx, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Vz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      if(Pml->getApplypml(0)){
         memcpy(this->checkpoints+pos, Pml->Sxx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(1)){
         memcpy(this->checkpoints+pos, Pml->Sxx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(4)){
         memcpy(this->checkpoints+pos, Pml->Szz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;
      }

      if(Pml->getApplypml(5)){
         memcpy(this->checkpoints+pos, Pml->Szz_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzz_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzz_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxz_bottom, nx_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->write(Sxx, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->write(Szz, nz_pml*nx_pml);
         this->Fc->write(Sxz, nz_pml*nx_pml);
         this->Fc->write(Vx, nz_pml*nx_pml);
         this->Fc->write(Vz, nz_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->write(Pml->Sxx_left, nz_pml*lpml);
            this->Fc->write(Pml->Sxzx_left, nz_pml*lpml);
            this->Fc->write(Pml->Vxx_left, nz_pml*lpml);
            this->Fc->write(Pml->Vzx_left, nz_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->write(Pml->Sxx_right, nz_pml*lpml);
            this->Fc->write(Pml->Sxzx_right, nz_pml*lpml);
            this->Fc->write(Pml->Vxx_right, nz_pml*lpml);
            this->Fc->write(Pml->Vzx_right, nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->write(Pml->Szz_top, nx_pml*lpml);
            this->Fc->write(Pml->Sxzz_top, nx_pml*lpml);
            this->Fc->write(Pml->Vzz_top, nx_pml*lpml);
            this->Fc->write(Pml->Vxz_top, nx_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->write(Pml->Szz_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Sxzz_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Vzz_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Vxz_bottom, nx_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::writeCheck: Error writing checkpoints to file.");
      }else{
         rs_error("Revolve::writeCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::readCheck(std::shared_ptr<WavesElastic3D<T>> waves)
{

   size_t nx_pml, ny_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Syy = waves->getSyy();
   T *Szz = waves->getSzz();
   T *Syz = waves->getSyz();
   T *Sxz = waves->getSxz();
   T *Sxy = waves->getSxy();
   T *Vx = waves->getVx();
   T *Vy = waves->getVy();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic3D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::readCheck: checkpoint array is not allocated.");
      memcpy(Sxx, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Syy, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Szz, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Syz, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Sxz, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Sxy, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Vx, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Vy, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Vz, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      if(Pml->getApplypml(0)){
         memcpy(Pml->Sxx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Sxzx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Sxyx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vxx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vzx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vyx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(1)){
         memcpy(Pml->Sxx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Sxzx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Sxyx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vxx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vzx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vyx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(2)){
         memcpy(Pml->Syy_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Sxyy_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Syzy_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Vyy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(Pml->Vxy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(Pml->Vzy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(3)){
         memcpy(Pml->Syy_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Sxyy_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Syzy_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Vyy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(Pml->Vxy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(Pml->Vzy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(4)){
         memcpy(Pml->Szz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Sxzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Syzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vxz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vyz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(5)){
         memcpy(Pml->Szz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Sxzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Syzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vxz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vyz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      }

   }else{
      if(this->open){
         this->Fc->read(Sxx,nz_pml*ny_pml*nx_pml, pos*sizeof(T));
         this->Fc->read(Syy,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Szz,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Syz,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Sxz,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Sxy,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Vx,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Vy,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Vz,nz_pml*ny_pml*nx_pml);

         if(Pml->getApplypml(0)){
            this->Fc->read(Pml->Sxx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxzx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxyx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vyx_left,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->read(Pml->Sxx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxzx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxyx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vyx_right,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(2)){
            this->Fc->read(Pml->Syy_front,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Sxyy_front,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Syzy_front,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Vyy_front,nx_pml*nz_pml*lpml);
            this->Fc->read(Pml->Vxy_front,nx_pml*nz_pml*lpml);
            this->Fc->read(Pml->Vzy_front,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(3)){
            this->Fc->read(Pml->Syy_back,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Sxyy_back,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Syzy_back,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Vyy_back,nx_pml*nz_pml*lpml);
            this->Fc->read(Pml->Vxy_back,nx_pml*nz_pml*lpml);
            this->Fc->read(Pml->Vzy_back,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->read(Pml->Szz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxzz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Syzz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vyz_top,nx_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->read(Pml->Szz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Syzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vyz_bottom,nx_pml*ny_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::readCheck: Error reading checkpoints from file.");
      }else{
         rs_error("Revolve::readCheck: File is closed.");

      }
   }
}

   template<typename T>
void Revolve<T>::writeCheck(std::shared_ptr<WavesElastic3D<T>> waves)
{

   size_t nx_pml, ny_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Syy = waves->getSyy();
   T *Szz = waves->getSzz();
   T *Syz = waves->getSyz();
   T *Sxz = waves->getSxz();
   T *Sxy = waves->getSxy();
   T *Vx = waves->getVx();
   T *Vy = waves->getVy();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic3D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::writeCheck: checkpoint array is not allocated.");
      memcpy(this->checkpoints+pos, Sxx, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Syy, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Szz, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Syz, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Sxz, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Sxy, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Vx, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Vy, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Vz, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      if(Pml->getApplypml(0)){
         memcpy(this->checkpoints+pos, Pml->Sxx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxzx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxyx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(1)){
         memcpy(this->checkpoints+pos, Pml->Sxx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxzx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxyx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(2)){
         memcpy(this->checkpoints+pos, Pml->Syy_front, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxyy_front, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Syzy_front, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyy_front, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxy_front, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzy_front, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(3)){
         memcpy(this->checkpoints+pos, Pml->Syy_back, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxyy_back, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Syzy_back, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyy_back, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxy_back, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzy_back, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(4)){
         memcpy(this->checkpoints+pos, Pml->Szz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxzz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Syzz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(5)){
         memcpy(this->checkpoints+pos, Pml->Szz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Syzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->write(Sxx,nz_pml*ny_pml*nx_pml, pos*sizeof(T));
         this->Fc->write(Syy,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Szz,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Syz,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Sxz,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Sxy,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Vx,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Vy,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Vz,nz_pml*ny_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->write(Pml->Sxx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxzx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxyx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vyx_left,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->write(Pml->Sxx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxzx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxyx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vyx_right,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(2)){
            this->Fc->write(Pml->Syy_front,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Sxyy_front,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Syzy_front,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Vyy_front,nx_pml*nz_pml*lpml);
            this->Fc->write(Pml->Vxy_front,nx_pml*nz_pml*lpml);
            this->Fc->write(Pml->Vzy_front,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(3)){
            this->Fc->write(Pml->Syy_back,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Sxyy_back,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Syzy_back,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Vyy_back,nx_pml*nz_pml*lpml);
            this->Fc->write(Pml->Vxy_back,nx_pml*nz_pml*lpml);
            this->Fc->write(Pml->Vzy_back,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->write(Pml->Szz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxzz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Syzz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vyz_top,nx_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->write(Pml->Szz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Syzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vyz_bottom,nx_pml*ny_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::writeCheck: Error writing checkpoints to file.");
      }else{
         rs_error("Revolve::writeCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::readCheck(std::shared_ptr<WavesElastic2D_DS<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Ux1 = waves->getUx1();
   T *Ux2 = waves->getUx2();
   T *Uz1 = waves->getUz1();
   T *Uz2 = waves->getUz2();
   std::shared_ptr<PmlElastic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::readCheck: checkpoint array is not allocated.");;
      memcpy(Ux1, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Ux2, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Uz1, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Uz2, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Pml->Sxx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(Pml->Sxx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(Pml->Sxzx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(Pml->Sxzx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(Pml->Szz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(Pml->Szz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(Pml->Sxzz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(Pml->Sxzz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(Pml->Vxx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(Pml->Vxx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(Pml->Vzx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(Pml->Vzx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(Pml->Vzz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(Pml->Vzz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(Pml->Vxz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(Pml->Vxz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
   }else{
      if(this->open){
         this->Fc->read(Ux1, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->read(Ux2, nz_pml*nx_pml);
         this->Fc->read(Uz1, nz_pml*nx_pml);
         this->Fc->read(Uz2, nz_pml*nx_pml);
         this->Fc->read(Pml->Sxx_left, nz_pml*lpml);
         this->Fc->read(Pml->Sxx_right, nz_pml*lpml);
         this->Fc->read(Pml->Sxzx_left, nz_pml*lpml);
         this->Fc->read(Pml->Sxzx_right, nz_pml*lpml);
         this->Fc->read(Pml->Szz_top, nx_pml*lpml);
         this->Fc->read(Pml->Szz_bottom, nx_pml*lpml);
         this->Fc->read(Pml->Sxzz_top, nx_pml*lpml);
         this->Fc->read(Pml->Sxzz_bottom, nx_pml*lpml);
         this->Fc->read(Pml->Vxx_left, nz_pml*lpml);
         this->Fc->read(Pml->Vxx_right, nz_pml*lpml);
         this->Fc->read(Pml->Vzx_left, nz_pml*lpml);
         this->Fc->read(Pml->Vzx_right, nz_pml*lpml);
         this->Fc->read(Pml->Vzz_top, nx_pml*lpml);
         this->Fc->read(Pml->Vzz_bottom, nx_pml*lpml);
         this->Fc->read(Pml->Vxz_top, nx_pml*lpml);
         this->Fc->read(Pml->Vxz_bottom, nx_pml*lpml);
         if(Fc->getFail()) rs_error("Revolve::readCheck: Error reading checkpoints from file.");
      }else{
         rs_error("Revolve::readCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::writeCheck(std::shared_ptr<WavesElastic2D_DS<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Ux1 = waves->getUx1();
   T *Ux2 = waves->getUx2();
   T *Uz1 = waves->getUz1();
   T *Uz2 = waves->getUz2();
   std::shared_ptr<PmlElastic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::writeCheck: checkpoint array is not allocated.");;

      memcpy(this->checkpoints+pos, Ux1, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Ux2, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Uz1, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Uz2, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Pml->Sxx_left, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Sxx_right, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Sxzx_left, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Sxzx_right, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Szz_top, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Szz_bottom, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Sxzz_top, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Sxzz_bottom, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Vxx_left, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Vxx_right, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Vzx_left, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Vzx_right, nz_pml*lpml*sizeof(T));
      pos += nz_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Vzz_top, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Vzz_bottom, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Vxz_top, nx_pml*lpml*sizeof(T));
      pos += nx_pml*lpml;

      memcpy(this->checkpoints+pos, Pml->Vxz_bottom, nx_pml*lpml*sizeof(T));
   }else{
      if(this->open){
         this->Fc->write(Ux1, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->write(Ux2, nz_pml*nx_pml);
         this->Fc->write(Uz1, nz_pml*nx_pml);
         this->Fc->write(Uz2, nz_pml*nx_pml);
         this->Fc->write(Pml->Sxx_left, nz_pml*lpml);
         this->Fc->write(Pml->Sxx_right, nz_pml*lpml);
         this->Fc->write(Pml->Sxzx_left, nz_pml*lpml);
         this->Fc->write(Pml->Sxzx_right, nz_pml*lpml);
         this->Fc->write(Pml->Szz_top, nx_pml*lpml);
         this->Fc->write(Pml->Szz_bottom, nx_pml*lpml);
         this->Fc->write(Pml->Sxzz_top, nx_pml*lpml);
         this->Fc->write(Pml->Sxzz_bottom, nx_pml*lpml);
         this->Fc->write(Pml->Vxx_left, nz_pml*lpml);
         this->Fc->write(Pml->Vxx_right, nz_pml*lpml);
         this->Fc->write(Pml->Vzx_left, nz_pml*lpml);
         this->Fc->write(Pml->Vzx_right, nz_pml*lpml);
         this->Fc->write(Pml->Vzz_top, nx_pml*lpml);
         this->Fc->write(Pml->Vzz_bottom, nx_pml*lpml);
         this->Fc->write(Pml->Vxz_top, nx_pml*lpml);
         this->Fc->write(Pml->Vxz_bottom, nx_pml*lpml);
         if(Fc->getFail()) rs_error("Revolve::writeCheck: Error writing checkpoints to file.");
      }else{
         rs_error("Revolve::writeCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::readCheck(std::shared_ptr<WavesElastic3D_DS<T>> waves)
{

   size_t nx_pml, ny_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Ux1 = waves->getUx1();
   T *Ux2 = waves->getUx2();
   T *Uy1 = waves->getUy1();
   T *Uy2 = waves->getUy2();
   T *Uz1 = waves->getUz1();
   T *Uz2 = waves->getUz2();
   std::shared_ptr<PmlElastic3D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::readCheck: checkpoint array is not allocated.");
      memcpy(Ux1, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Ux2, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Uy1, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Uy2, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Uz1, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Uz2, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Pml->Sxx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Sxx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Sxzx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Sxzx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Sxyx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Sxyx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Vxx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Vxx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Vzx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Vzx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Vyx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Vyx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(Pml->Szz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Szz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Sxzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Sxzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Syzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Syzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Vzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Vzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Vxz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Vxz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Vyz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Vyz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(Pml->Syy_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(Pml->Syy_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(Pml->Sxyy_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(Pml->Sxyy_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(Pml->Syzy_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(Pml->Syzy_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(Pml->Vyy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(Pml->Vyy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(Pml->Vxy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(Pml->Vxy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(Pml->Vzy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(Pml->Vzy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
   }else{
      if(this->open){
         this->Fc->read(Ux1,nz_pml*ny_pml*nx_pml, pos*sizeof(T));
         this->Fc->read(Ux2,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Uy1,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Uy2,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Uz1,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Uz2,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Pml->Sxx_left,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Sxx_right,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Sxzx_left,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Sxzx_right,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Sxyx_left,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Sxyx_right,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vxx_left,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vxx_right,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vzx_left,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vzx_right,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vyx_left,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vyx_right,nz_pml*ny_pml*lpml);
         this->Fc->read(Pml->Syy_front,nz_pml*nx_pml*lpml);
         this->Fc->read(Pml->Syy_back,nz_pml*nx_pml*lpml);
         this->Fc->read(Pml->Sxyy_front,nz_pml*nx_pml*lpml);
         this->Fc->read(Pml->Sxyy_back,nz_pml*nx_pml*lpml);
         this->Fc->read(Pml->Syzy_front,nz_pml*nx_pml*lpml);
         this->Fc->read(Pml->Syzy_back,nz_pml*nx_pml*lpml);
         this->Fc->read(Pml->Vyy_front,nx_pml*nz_pml*lpml);
         this->Fc->read(Pml->Vyy_back,nx_pml*nz_pml*lpml);
         this->Fc->read(Pml->Vxy_front,nx_pml*nz_pml*lpml);
         this->Fc->read(Pml->Vxy_back,nx_pml*nz_pml*lpml);
         this->Fc->read(Pml->Vzy_front,nx_pml*nz_pml*lpml);
         this->Fc->read(Pml->Vzy_back,nx_pml*nz_pml*lpml);
         this->Fc->read(Pml->Szz_top,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Szz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Sxzz_top,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Sxzz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Syzz_top,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Syzz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vzz_top,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vzz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vxz_top,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vxz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vyz_top,nx_pml*ny_pml*lpml);
         this->Fc->read(Pml->Vyz_bottom,nx_pml*ny_pml*lpml);
         if(Fc->getFail()) rs_error("Revolve::readCheck: Error reading checkpoints from file.");
      }else{
         rs_error("Revolve::readCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::writeCheck(std::shared_ptr<WavesElastic3D_DS<T>> waves)
{

   size_t nx_pml, ny_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Ux1 = waves->getUx1();
   T *Ux2 = waves->getUx2();
   T *Uy1 = waves->getUy1();
   T *Uy2 = waves->getUy2();
   T *Uz1 = waves->getUz1();
   T *Uz2 = waves->getUz2();
   std::shared_ptr<PmlElastic3D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::writeCheck: checkpoint array is not allocated.");
      memcpy(this->checkpoints+pos, Ux1, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Ux2, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Uy1, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Uy2, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Uz1, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Uz2, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Pml->Sxx_left, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Sxx_right, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Sxzx_left, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Sxzx_right, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Sxyx_left, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Sxyx_right, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vxx_left, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vxx_right, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vzx_left, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vzx_right, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vyx_left, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vyx_right, nz_pml*ny_pml*lpml*sizeof(T));
      pos += nz_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Szz_top, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Szz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Sxzz_top, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Sxzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Syzz_top, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Syzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vzz_top, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vxz_top, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vxz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vyz_top, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vyz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
      pos += nx_pml*ny_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Syy_front, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Syy_back, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Sxyy_front, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Sxyy_back, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Syzy_front, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Syzy_back, nz_pml*nx_pml*lpml*sizeof(T));
      pos += nz_pml*nx_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vyy_front, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vyy_back, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vxy_front, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vxy_back, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vzy_front, nx_pml*nz_pml*lpml*sizeof(T));
      pos += nx_pml*nz_pml*lpml;
      memcpy(this->checkpoints+pos, Pml->Vzy_back, nx_pml*nz_pml*lpml*sizeof(T));
   }else{
      if(this->open){
         this->Fc->write(Ux1,nz_pml*ny_pml*nx_pml, pos*sizeof(T));
         this->Fc->write(Ux2,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Uy1,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Uy2,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Uz1,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Uz2,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Pml->Sxx_left,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Sxx_right,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Sxzx_left,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Sxzx_right,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Sxyx_left,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Sxyx_right,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vxx_left,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vxx_right,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vzx_left,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vzx_right,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vyx_left,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vyx_right,nz_pml*ny_pml*lpml);
         this->Fc->write(Pml->Syy_front,nz_pml*nx_pml*lpml);
         this->Fc->write(Pml->Syy_back,nz_pml*nx_pml*lpml);
         this->Fc->write(Pml->Sxyy_front,nz_pml*nx_pml*lpml);
         this->Fc->write(Pml->Sxyy_back,nz_pml*nx_pml*lpml);
         this->Fc->write(Pml->Syzy_front,nz_pml*nx_pml*lpml);
         this->Fc->write(Pml->Syzy_back,nz_pml*nx_pml*lpml);
         this->Fc->write(Pml->Vyy_front,nx_pml*nz_pml*lpml);
         this->Fc->write(Pml->Vyy_back,nx_pml*nz_pml*lpml);
         this->Fc->write(Pml->Vxy_front,nx_pml*nz_pml*lpml);
         this->Fc->write(Pml->Vxy_back,nx_pml*nz_pml*lpml);
         this->Fc->write(Pml->Vzy_front,nx_pml*nz_pml*lpml);
         this->Fc->write(Pml->Vzy_back,nx_pml*nz_pml*lpml);
         this->Fc->write(Pml->Szz_top,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Szz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Sxzz_top,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Sxzz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Syzz_top,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Syzz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vzz_top,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vzz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vxz_top,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vxz_bottom,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vyz_top,nx_pml*ny_pml*lpml);
         this->Fc->write(Pml->Vyz_bottom,nx_pml*ny_pml*lpml);
         if(Fc->getFail()) rs_error("Revolve::writeCheck: Error writing checkpoints to file.");
      }else{
         rs_error("Revolve::writeCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::readCheck(std::shared_ptr<WavesViscoelastic2D<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Szz = waves->getSzz();
   T *Sxz = waves->getSxz();
   T *Mxx = waves->getMxx();
   T *Mzz = waves->getMzz();
   T *Mxz = waves->getMxz();
   T *Vx = waves->getVx();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::readCheck: checkpoint array is not allocated.");;
      memcpy(Sxx, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Szz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Sxz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Mxx, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Mzz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Mxz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Vx, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Vz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      if(Pml->getApplypml(0)){
         memcpy(Pml->Sxx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Sxzx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vxx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vzx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(1)){
         memcpy(Pml->Sxx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Sxzx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vxx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vzx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(4)){
         memcpy(Pml->Szz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Sxzz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vzz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vxz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;
      }

      if(Pml->getApplypml(5)){
         memcpy(Pml->Szz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Sxzz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vzz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vxz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->read(Sxx, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->read(Szz, nz_pml*nx_pml);
         this->Fc->read(Sxz, nz_pml*nx_pml);
         this->Fc->read(Mxx, nz_pml*nx_pml);
         this->Fc->read(Mzz, nz_pml*nx_pml);
         this->Fc->read(Mxz, nz_pml*nx_pml);
         this->Fc->read(Vx, nz_pml*nx_pml);
         this->Fc->read(Vz, nz_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->read(Pml->Sxx_left, nz_pml*lpml);
            this->Fc->read(Pml->Sxzx_left, nz_pml*lpml);
            this->Fc->read(Pml->Vxx_left, nz_pml*lpml);
            this->Fc->read(Pml->Vzx_left, nz_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->read(Pml->Sxx_right, nz_pml*lpml);
            this->Fc->read(Pml->Sxzx_right, nz_pml*lpml);
            this->Fc->read(Pml->Vxx_right, nz_pml*lpml);
            this->Fc->read(Pml->Vzx_right, nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->read(Pml->Szz_top, nx_pml*lpml);
            this->Fc->read(Pml->Sxzz_top, nx_pml*lpml);
            this->Fc->read(Pml->Vzz_top, nx_pml*lpml);
            this->Fc->read(Pml->Vxz_top, nx_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->read(Pml->Szz_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Sxzz_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Vzz_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Vxz_bottom, nx_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::readCheck: Error reading checkpoints from file.");
      }else{
         rs_error("Revolve::readCheck: File is closed.");

      }
   }
}

   template<typename T>
void Revolve<T>::writeCheck(std::shared_ptr<WavesViscoelastic2D<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Szz = waves->getSzz();
   T *Sxz = waves->getSxz();
   T *Mxx = waves->getMxx();
   T *Mzz = waves->getMzz();
   T *Mxz = waves->getMxz();
   T *Vx = waves->getVx();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::writeCheck: checkpoint array is not allocated.");;
      memcpy(this->checkpoints+pos, Sxx, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Szz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Sxz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Mxx, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Mzz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Mxz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Vx, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Vz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      if(Pml->getApplypml(0)){
         memcpy(this->checkpoints+pos, Pml->Sxx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(1)){
         memcpy(this->checkpoints+pos, Pml->Sxx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(4)){
         memcpy(this->checkpoints+pos, Pml->Szz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;
      }

      if(Pml->getApplypml(5)){
         memcpy(this->checkpoints+pos, Pml->Szz_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzz_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzz_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxz_bottom, nx_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->write(Sxx, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->write(Szz, nz_pml*nx_pml);
         this->Fc->write(Sxz, nz_pml*nx_pml);
         this->Fc->write(Mxx, nz_pml*nx_pml);
         this->Fc->write(Mzz, nz_pml*nx_pml);
         this->Fc->write(Mxz, nz_pml*nx_pml);
         this->Fc->write(Vx, nz_pml*nx_pml);
         this->Fc->write(Vz, nz_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->write(Pml->Sxx_left, nz_pml*lpml);
            this->Fc->write(Pml->Sxzx_left, nz_pml*lpml);
            this->Fc->write(Pml->Vxx_left, nz_pml*lpml);
            this->Fc->write(Pml->Vzx_left, nz_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->write(Pml->Sxx_right, nz_pml*lpml);
            this->Fc->write(Pml->Sxzx_right, nz_pml*lpml);
            this->Fc->write(Pml->Vxx_right, nz_pml*lpml);
            this->Fc->write(Pml->Vzx_right, nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->write(Pml->Szz_top, nx_pml*lpml);
            this->Fc->write(Pml->Sxzz_top, nx_pml*lpml);
            this->Fc->write(Pml->Vzz_top, nx_pml*lpml);
            this->Fc->write(Pml->Vxz_top, nx_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->write(Pml->Szz_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Sxzz_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Vzz_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Vxz_bottom, nx_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::writeCheck: Error writing checkpoints to file.");
      }else{
         rs_error("Revolve::writeCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::readCheck(std::shared_ptr<WavesVti2D<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Szz = waves->getSzz();
   T *Sxz = waves->getSxz();
   T *Vx = waves->getVx();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::readCheck: checkpoint array is not allocated.");;
      memcpy(Sxx, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Szz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Sxz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Vx, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(Vz, this->checkpoints+pos, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      if(Pml->getApplypml(0)){
         memcpy(Pml->Sxx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Sxzx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vxx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vzx_left, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(1)){
         memcpy(Pml->Sxx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Sxzx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vxx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(Pml->Vzx_right, this->checkpoints+pos, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(4)){
         memcpy(Pml->Szz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Sxzz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vzz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vxz_top, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;
      }

      if(Pml->getApplypml(5)){
         memcpy(Pml->Szz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Sxzz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vzz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(Pml->Vxz_bottom, this->checkpoints+pos, nx_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->read(Sxx, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->read(Szz, nz_pml*nx_pml);
         this->Fc->read(Sxz, nz_pml*nx_pml);
         this->Fc->read(Vx, nz_pml*nx_pml);
         this->Fc->read(Vz, nz_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->read(Pml->Sxx_left, nz_pml*lpml);
            this->Fc->read(Pml->Sxzx_left, nz_pml*lpml);
            this->Fc->read(Pml->Vxx_left, nz_pml*lpml);
            this->Fc->read(Pml->Vzx_left, nz_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->read(Pml->Sxx_right, nz_pml*lpml);
            this->Fc->read(Pml->Sxzx_right, nz_pml*lpml);
            this->Fc->read(Pml->Vxx_right, nz_pml*lpml);
            this->Fc->read(Pml->Vzx_right, nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->read(Pml->Szz_top, nx_pml*lpml);
            this->Fc->read(Pml->Sxzz_top, nx_pml*lpml);
            this->Fc->read(Pml->Vzz_top, nx_pml*lpml);
            this->Fc->read(Pml->Vxz_top, nx_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->read(Pml->Szz_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Sxzz_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Vzz_bottom, nx_pml*lpml);
            this->Fc->read(Pml->Vxz_bottom, nx_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::readCheck: Error reading checkpoints from file.");
      }else{
         rs_error("Revolve::readCheck: File is closed.");

      }
   }
}

   template<typename T>
void Revolve<T>::writeCheck(std::shared_ptr<WavesVti2D<T>> waves)
{

   size_t nx_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Szz = waves->getSzz();
   T *Sxz = waves->getSxz();
   T *Vx = waves->getVx();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic2D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::writeCheck: checkpoint array is not allocated.");;
      memcpy(this->checkpoints+pos, Sxx, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Szz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Sxz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Vx, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      memcpy(this->checkpoints+pos, Vz, nz_pml*nx_pml*sizeof(T));
      pos += nz_pml*nx_pml;

      if(Pml->getApplypml(0)){
         memcpy(this->checkpoints+pos, Pml->Sxx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzx_left, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(1)){
         memcpy(this->checkpoints+pos, Pml->Sxx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzx_right, nz_pml*lpml*sizeof(T));
         pos += nz_pml*lpml;
      }

      if(Pml->getApplypml(4)){
         memcpy(this->checkpoints+pos, Pml->Szz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxz_top, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;
      }

      if(Pml->getApplypml(5)){
         memcpy(this->checkpoints+pos, Pml->Szz_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Sxzz_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vzz_bottom, nx_pml*lpml*sizeof(T));
         pos += nx_pml*lpml;

         memcpy(this->checkpoints+pos, Pml->Vxz_bottom, nx_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->write(Sxx, nz_pml*nx_pml, pos*sizeof(T));
         this->Fc->write(Szz, nz_pml*nx_pml);
         this->Fc->write(Sxz, nz_pml*nx_pml);
         this->Fc->write(Vx, nz_pml*nx_pml);
         this->Fc->write(Vz, nz_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->write(Pml->Sxx_left, nz_pml*lpml);
            this->Fc->write(Pml->Sxzx_left, nz_pml*lpml);
            this->Fc->write(Pml->Vxx_left, nz_pml*lpml);
            this->Fc->write(Pml->Vzx_left, nz_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->write(Pml->Sxx_right, nz_pml*lpml);
            this->Fc->write(Pml->Sxzx_right, nz_pml*lpml);
            this->Fc->write(Pml->Vxx_right, nz_pml*lpml);
            this->Fc->write(Pml->Vzx_right, nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->write(Pml->Szz_top, nx_pml*lpml);
            this->Fc->write(Pml->Sxzz_top, nx_pml*lpml);
            this->Fc->write(Pml->Vzz_top, nx_pml*lpml);
            this->Fc->write(Pml->Vxz_top, nx_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->write(Pml->Szz_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Sxzz_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Vzz_bottom, nx_pml*lpml);
            this->Fc->write(Pml->Vxz_bottom, nx_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::writeCheck: Error writing checkpoints to file.");
      }else{
         rs_error("Revolve::writeCheck: File is closed.");

      }
   }
}

template<typename T>
void Revolve<T>::readCheck(std::shared_ptr<WavesOrtho3D<T>> waves)
{

   size_t nx_pml, ny_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Syy = waves->getSyy();
   T *Szz = waves->getSzz();
   T *Syz = waves->getSyz();
   T *Sxz = waves->getSxz();
   T *Sxy = waves->getSxy();
   T *Vx = waves->getVx();
   T *Vy = waves->getVy();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic3D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::readCheck: checkpoint array is not allocated.");
      memcpy(Sxx, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Syy, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Szz, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Syz, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Sxz, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Sxy, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Vx, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Vy, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(Vz, this->checkpoints+pos, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      if(Pml->getApplypml(0)){
         memcpy(Pml->Sxx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Sxzx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Sxyx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vxx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vzx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vyx_left, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(1)){
         memcpy(Pml->Sxx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Sxzx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Sxyx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vxx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vzx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(Pml->Vyx_right, this->checkpoints+pos, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(2)){
         memcpy(Pml->Syy_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Sxyy_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Syzy_front, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Vyy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(Pml->Vxy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(Pml->Vzy_front, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(3)){
         memcpy(Pml->Syy_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Sxyy_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Syzy_back, this->checkpoints+pos, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(Pml->Vyy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(Pml->Vxy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(Pml->Vzy_back, this->checkpoints+pos, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(4)){
         memcpy(Pml->Szz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Sxzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Syzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vzz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vxz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vyz_top, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(5)){
         memcpy(Pml->Szz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Sxzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Syzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vzz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vxz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(Pml->Vyz_bottom, this->checkpoints+pos, nx_pml*ny_pml*lpml*sizeof(T));
      }

   }else{
      if(this->open){
         this->Fc->read(Sxx,nz_pml*ny_pml*nx_pml, pos*sizeof(T));
         this->Fc->read(Syy,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Szz,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Syz,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Sxz,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Sxy,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Vx,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Vy,nz_pml*ny_pml*nx_pml);
         this->Fc->read(Vz,nz_pml*ny_pml*nx_pml);

         if(Pml->getApplypml(0)){
            this->Fc->read(Pml->Sxx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxzx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxyx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzx_left,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vyx_left,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->read(Pml->Sxx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxzx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxyx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzx_right,nz_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vyx_right,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(2)){
            this->Fc->read(Pml->Syy_front,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Sxyy_front,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Syzy_front,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Vyy_front,nx_pml*nz_pml*lpml);
            this->Fc->read(Pml->Vxy_front,nx_pml*nz_pml*lpml);
            this->Fc->read(Pml->Vzy_front,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(3)){
            this->Fc->read(Pml->Syy_back,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Sxyy_back,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Syzy_back,nz_pml*nx_pml*lpml);
            this->Fc->read(Pml->Vyy_back,nx_pml*nz_pml*lpml);
            this->Fc->read(Pml->Vxy_back,nx_pml*nz_pml*lpml);
            this->Fc->read(Pml->Vzy_back,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->read(Pml->Szz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxzz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Syzz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxz_top,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vyz_top,nx_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->read(Pml->Szz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Sxzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Syzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vxz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->read(Pml->Vyz_bottom,nx_pml*ny_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::readCheck: Error reading checkpoints from file.");
      }else{
         rs_error("Revolve::readCheck: File is closed.");

      }
   }
}

   template<typename T>
void Revolve<T>::writeCheck(std::shared_ptr<WavesOrtho3D<T>> waves)
{

   size_t nx_pml, ny_pml, nz_pml;
   size_t lpml;
   nx_pml=waves->getNx_pml();
   ny_pml=waves->getNy_pml();
   nz_pml=waves->getNz_pml();
   lpml = waves->getLpml();
   T *Sxx = waves->getSxx();
   T *Syy = waves->getSyy();
   T *Szz = waves->getSzz();
   T *Syz = waves->getSyz();
   T *Sxz = waves->getSxz();
   T *Sxy = waves->getSxy();
   T *Vx = waves->getVx();
   T *Vy = waves->getVy();
   T *Vz = waves->getVz();
   std::shared_ptr<PmlElastic3D<T>> Pml = waves->getPml();

   off_t pos = this->checksize*this->check; 
   if(this->incore){
      if(!this->allocated) rs_error("Revolve::writeCheck: checkpoint array is not allocated.");
      memcpy(this->checkpoints+pos, Sxx, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Syy, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Szz, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Syz, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Sxz, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Sxy, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Vx, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Vy, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      memcpy(this->checkpoints+pos, Vz, nz_pml*ny_pml*nx_pml*sizeof(T));
      pos += nz_pml*ny_pml*nx_pml;
      if(Pml->getApplypml(0)){
         memcpy(this->checkpoints+pos, Pml->Sxx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxzx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxyx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyx_left, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(1)){
         memcpy(this->checkpoints+pos, Pml->Sxx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxzx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxyx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyx_right, nz_pml*ny_pml*lpml*sizeof(T));
         pos += nz_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(2)){
         memcpy(this->checkpoints+pos, Pml->Syy_front, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxyy_front, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Syzy_front, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyy_front, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxy_front, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzy_front, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(3)){
         memcpy(this->checkpoints+pos, Pml->Syy_back, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxyy_back, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Syzy_back, nz_pml*nx_pml*lpml*sizeof(T));
         pos += nz_pml*nx_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyy_back, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxy_back, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzy_back, nx_pml*nz_pml*lpml*sizeof(T));
         pos += nx_pml*nz_pml*lpml;
      }
      if(Pml->getApplypml(4)){
         memcpy(this->checkpoints+pos, Pml->Szz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxzz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Syzz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyz_top, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
      }
      if(Pml->getApplypml(5)){
         memcpy(this->checkpoints+pos, Pml->Szz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Sxzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Syzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vzz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vxz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
         pos += nx_pml*ny_pml*lpml;
         memcpy(this->checkpoints+pos, Pml->Vyz_bottom, nx_pml*ny_pml*lpml*sizeof(T));
      }
   }else{
      if(this->open){
         this->Fc->write(Sxx,nz_pml*ny_pml*nx_pml, pos*sizeof(T));
         this->Fc->write(Syy,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Szz,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Syz,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Sxz,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Sxy,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Vx,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Vy,nz_pml*ny_pml*nx_pml);
         this->Fc->write(Vz,nz_pml*ny_pml*nx_pml);
         if(Pml->getApplypml(0)){
            this->Fc->write(Pml->Sxx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxzx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxyx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzx_left,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vyx_left,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(1)){
            this->Fc->write(Pml->Sxx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxzx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxyx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzx_right,nz_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vyx_right,nz_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(2)){
            this->Fc->write(Pml->Syy_front,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Sxyy_front,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Syzy_front,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Vyy_front,nx_pml*nz_pml*lpml);
            this->Fc->write(Pml->Vxy_front,nx_pml*nz_pml*lpml);
            this->Fc->write(Pml->Vzy_front,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(3)){
            this->Fc->write(Pml->Syy_back,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Sxyy_back,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Syzy_back,nz_pml*nx_pml*lpml);
            this->Fc->write(Pml->Vyy_back,nx_pml*nz_pml*lpml);
            this->Fc->write(Pml->Vxy_back,nx_pml*nz_pml*lpml);
            this->Fc->write(Pml->Vzy_back,nx_pml*nz_pml*lpml);
         }
         if(Pml->getApplypml(4)){
            this->Fc->write(Pml->Szz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxzz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Syzz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxz_top,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vyz_top,nx_pml*ny_pml*lpml);
         }
         if(Pml->getApplypml(5)){
            this->Fc->write(Pml->Szz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Sxzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Syzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vzz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vxz_bottom,nx_pml*ny_pml*lpml);
            this->Fc->write(Pml->Vyz_bottom,nx_pml*ny_pml*lpml);
         }
         if(Fc->getFail()) rs_error("Revolve::writeCheck: Error writing checkpoints to file.");
      }else{
         rs_error("Revolve::writeCheck: File is closed.");

      }
   }
}

// destructor
template<typename T>
Revolve<T>::~Revolve(){
   if(this->allocated){
      free(checkpoints);
   }
   if(this->open){
      Fc->close();
   }
}

// =============== INITIALIZING TEMPLATE STRUCTS =============== //
template class Revolve<float>;
template class Revolve<double>;

}


