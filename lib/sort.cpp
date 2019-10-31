#include "sort.h"

namespace rockseis {
// constructor
template<typename T>
Sort<T>::Sort()
{
    ngathers = 0;
    ntraces = 0;
    sortkey = SOURCE;
    kmapfile = "kmap.rss";
    smapfile = "smap.rss";
    keymap=(key *) calloc(1,1);
    sortmap=(size_t *) calloc(1,1);
    reciprocity = false;
    sortset = false;
}

template<typename T>
Sort<T>::~Sort() {
    free(keymap);
    free(sortmap);
}

//Sort with x as primary key and y as secondary key 
int sort_sr_0(void const *a, void const *b)
{
	position_t *pa, *pb;

	pa = (position_t *) a;
	pb = (position_t *) b;

	if( pa->z < pb->z) return -1;
	if( pa->z > pb->z) return 1;

	if( pa->y < pb->y) return -1;
	if( pa->y > pb->y) return 1;

	if( pa->x < pb->x) return -1;
	if( pa->x > pb->x) return 1;

	if( pa->foff < pb->foff) return -1;
	if( pa->foff > pb->foff) return 1;

	if( pa->offz < pb->offz) return -1;
	if( pa->offz > pb->offz) return 1;

	if( pa->offy < pb->offy) return -1;
	if( pa->offy > pb->offy) return 1;

	if( pa->offx < pb->offx) return -1;
	if( pa->offx > pb->offx) return 1;

	return 0;
}

template<typename T>
bool Sort<T>::createSort(std::string filename, rs_key _sortkey, T dx, T dy)
{
    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    if(strcmp(filename.c_str(), "stdin")){
        status = Fdata->input(filename);
    }else{
        status = Fdata->input();
    }
    if(status == FILE_ERR) rs_error("Sort::createSort: Error reading input file: ", filename);
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA2D && datatype != DATA3D) rs_error("Sort::createSort: Only DATA2D and DATA3D types are supported.");

    sortkey = _sortkey; 
    datafile = filename;
    size_t n1 = Fdata->getN(1);
    size_t n2 = Fdata->getN(2);

    int nh = Fdata->getNheader();
    int dp = Fdata->getData_format();
    int hp = Fdata->getHeader_format();
    float fbuffer[NHEAD3D];
	double dbuffer[NHEAD3D];
	float ftmp;
	double dtmp;
    position_t *positions, *ensembles;

	positions = (position_t *) calloc(n2,sizeof(position_t));
	ensembles = (position_t *) calloc(n2,sizeof(position_t));

    size_t nensembles = 0;
    for (size_t i2=0; i2 < n2; i2++){
        //Read coordinates
        Fdata->seekg(i2*(n1*dp + nh*hp) + Fdata->getStartofdata());
        switch (hp){
            case sizeof(float):
                Fdata->read(fbuffer, nh);
                switch (sortkey){
                    case SOURCE:
                        if(datatype == DATA3D) {
                            positions[i2].x = fbuffer[0];
                            positions[i2].offx = fbuffer[3] - fbuffer[0];
                            positions[i2].y = fbuffer[1];
                            positions[i2].offy = fbuffer[4] - fbuffer[1];
                            positions[i2].z = fbuffer[2];
                            positions[i2].offz = fbuffer[5] - fbuffer[2];
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }else{
                            positions[i2].x = fbuffer[0];
                            positions[i2].offx = fbuffer[2] - fbuffer[0];
                            positions[i2].y = 0.0;
                            positions[i2].offy = 0.0;
                            positions[i2].z = fbuffer[1];
                            positions[i2].offz = fbuffer[3] - fbuffer[1];
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }
                        break;
                    case RECEIVER:
                        if(datatype == DATA3D) {
                            positions[i2].x = fbuffer[3];
                            positions[i2].offx = fbuffer[0] - fbuffer[3];
                            positions[i2].y = fbuffer[4];
                            positions[i2].offy = fbuffer[1] - fbuffer[4];
                            positions[i2].z = fbuffer[5];
                            positions[i2].offz = fbuffer[2] - fbuffer[5];
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }else{
                            positions[i2].x = fbuffer[2];
                            positions[i2].offx = fbuffer[0] - fbuffer[2];
                            positions[i2].y = 0.0;
                            positions[i2].offy = 0.0;
                            positions[i2].z = fbuffer[3];
                            positions[i2].offz = fbuffer[1] - fbuffer[3];
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }
                        break;
                    case CMP:
                        if(datatype == DATA3D) {
                            positions[i2].x = 0.5*(fbuffer[0]+fbuffer[3]);
                            positions[i2].offx = (fbuffer[3]-fbuffer[0]);
                            positions[i2].y = 0.5*(fbuffer[1]+fbuffer[4]);
                            positions[i2].offy = (fbuffer[4]-fbuffer[1]);
                            positions[i2].z = 0.0;
                            positions[i2].offz = 0.0;
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }else{
                            positions[i2].x = 0.5*(fbuffer[0]+fbuffer[2]);
                            positions[i2].offx = (fbuffer[2]-fbuffer[1]);
                            positions[i2].y = 0.0;
                            positions[i2].offy = 0.0;
                            positions[i2].z = 0.0;
                            positions[i2].offz = 0.0;
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }
                        // Snapping coordinates to bins
						ftmp=(positions[i2].x)/dx; // Compute index in irregular coordinate
						positions[i2].x=rintf(ftmp)*dx; // Compute nearest position in regular grid
						ftmp=(positions[i2].y)/dy; // Compute indey in irregular coordinate
						positions[i2].y=rintf(ftmp)*dy; // Compute nearest position in regular grid
						break;
                    default:
                        rs_error("Sort: sort key invalid.");
                        break;
                }
                break;
            case sizeof(double):
                Fdata->read(dbuffer, nh);
                switch (sortkey){
                    case SOURCE:
                        if(datatype == DATA3D) {
                            positions[i2].x = dbuffer[0];
                            positions[i2].offx = dbuffer[3] - dbuffer[0];
                            positions[i2].y = dbuffer[1];
                            positions[i2].offy = dbuffer[4] - dbuffer[1];
                            positions[i2].z = dbuffer[2];
                            positions[i2].offz = dbuffer[5] - dbuffer[2];
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }else{
                            positions[i2].x = dbuffer[0];
                            positions[i2].offx = dbuffer[2] - dbuffer[0];
                            positions[i2].y = 0.0;
                            positions[i2].offy = 0.0;
                            positions[i2].z = dbuffer[1];
                            positions[i2].offz = dbuffer[3] - dbuffer[1];
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }
                        break;
                    case RECEIVER:
                        if(datatype == DATA3D) {
                            positions[i2].x = dbuffer[3];
                            positions[i2].offx = dbuffer[0] - dbuffer[3];
                            positions[i2].y = dbuffer[4];
                            positions[i2].offy = dbuffer[1] - dbuffer[4];
                            positions[i2].z = dbuffer[5];
                            positions[i2].offz = dbuffer[2] - dbuffer[5];
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }else{
                            positions[i2].x = dbuffer[2];
                            positions[i2].offx = dbuffer[0] - dbuffer[2];
                            positions[i2].y = 0.0;
                            positions[i2].offy = 0.0;
                            positions[i2].z = dbuffer[3];
                            positions[i2].offz = dbuffer[1] - dbuffer[3];
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }
                        break;
                    case CMP:
                        if(datatype == DATA3D) {
                            positions[i2].x = 0.5*(dbuffer[0]+dbuffer[3]);
                            positions[i2].offx = (dbuffer[3]-dbuffer[0]);
                            positions[i2].y = 0.5*(dbuffer[1]+dbuffer[4]);
                            positions[i2].offy = (dbuffer[4]-dbuffer[1]);
                            positions[i2].z = 0.0;
                            positions[i2].offz = 0.0;
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }else{
                            positions[i2].x = 0.5*(dbuffer[0]+dbuffer[2]);
                            positions[i2].offx = (dbuffer[2]-dbuffer[1]);
                            positions[i2].y = 0.0;
                            positions[i2].offy = 0.0;
                            positions[i2].z = 0.0;
                            positions[i2].offz = 0.0;
                            positions[i2].foff = SQ(positions[i2].offx) + SQ(positions[i2].offy) + SQ(positions[i2].offz);
                        }
						// Snapping coordinates to bins
						dtmp=(positions[i2].x)/dx; // Compute index in irregular coordinate
						positions[i2].x=rintf(dtmp)*dx; // Compute nearest position in regular grid
						dtmp=(positions[i2].y)/dy; // Compute indey in irregular coordinate
						positions[i2].y=rintf(dtmp)*dy; // Compute nearest position in regular grid
						break;
                    default:
                        rs_error("Sort: sort key invalid.");
                        break;
                }
                break;
            default:
                rs_error("Sort::createSort: Header format not supported.");
                break;
        }
        positions[i2].ind = i2;
    }

	// Sort list
	qsort(positions, n2, sizeof(position_t), sort_sr_0);

    nensembles = 1;
    size_t ntraces = 1;
    double oldx, oldy, oldz;
    size_t oldind;

	/* Initialize the ensemble list and map */
	oldx = positions[0].x;
	oldy = positions[0].y;
	oldz = positions[0].z;
	ensembles[0].x=oldx;
	ensembles[0].y=oldy;
	ensembles[0].z=oldz;
	oldind = 0;

	// Making a map of the traces 
	size_t *map = (size_t *) calloc(2*n2,sizeof(size_t));
    Index I(2,n2);

	for (int i2=1; i2 < n2; i2++) {
		// Get coordinates
		dbuffer[0] = positions[i2].x;
		dbuffer[1] = positions[i2].y;
		dbuffer[2] = positions[i2].z;
		// If ensemble position is not in the list then add it to the list
		if((dbuffer[0] == oldx) && (dbuffer[1] == oldy) && (dbuffer[2] == oldz)){
			ntraces++;
		}else{
			ensembles[nensembles].x=dbuffer[0];
			ensembles[nensembles].y=dbuffer[1];
			ensembles[nensembles].z=dbuffer[2];
			map[I(0,nensembles-1)]=oldind;
			map[I(1,nensembles-1)]=ntraces;
			nensembles++;
			ntraces=1;
			oldx=dbuffer[0];
			oldy=dbuffer[1];
			oldz=dbuffer[2];
			oldind=i2;
		}
    }
    // Taking the last ensemble into consideration
	map[I(0,nensembles-1)]=oldind;
	map[I(1,nensembles-1)]=ntraces;

    // Create keymap and sortmap
    free(keymap); free(sortmap);
    keymap = (key *) calloc(nensembles, sizeof(key));
    sortmap = (size_t *) calloc(n2, sizeof(size_t));
    //Setting variables
    this->ngathers = nensembles;
    this->ntraces = n2;

	for (int i2=0; i2 < nensembles; i2++) {
        keymap[i2].i0 = map[I(0,i2)];
        keymap[i2].n = map[I(1,i2)];
        keymap[i2].status = NOT_STARTED;
    }

	for (int i2=0; i2 < n2; i2++) {
        sortmap[i2] = positions[i2].ind;
    }

    // Free allocated arrays
    free(positions);
    free(ensembles);
    free(map);

    sortset = true;
    return SORT_OK;
}

template<typename T>
std::shared_ptr<Data2D<T>> Sort<T>::get2DGather()
{
    if(this->ngathers == 0 || this->ntraces == 0) rs_error("No sort map created.");

    size_t i;
    for(i=0; i < this->ngathers; i++)
    {
        if(keymap[i].status == NOT_STARTED || keymap[i].status == FAILED)
        {
            break;
        }
    }
    if(i < this->ngathers) 
    {
        // Found a shot
        bool status;
        std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
        if(strcmp(this->datafile.c_str(), "stdin")){
            status = Fdata->input(this->datafile);
        }else{
            status = Fdata->input();
        }
        if(status == FILE_ERR) rs_error("Sort::get2DGather: Error reading input data file: ", this->datafile);
        rs_datatype datatype = Fdata->getType(); 
        if(datatype != DATA2D) rs_error("Sort::get2DGather: Datafile must be of type Data2D.");
        //Get gather size information
        size_t n1 = Fdata->getN(1);
        T d1 = Fdata->getD(1);
        T o1 = Fdata->getO(1);
        size_t n2 = this->keymap[i].n;

        //Create gather
	    std::shared_ptr<rockseis::Data2D<T>> gather;
        gather = std::make_shared<Data2D<T>>(n2,n1,d1,o1);
        Point2D<T> *scoords = (gather->getGeom())->getScoords();
        Point2D<T> *gcoords = (gather->getGeom())->getGcoords();
        T *data = gather->getData();
        size_t traceno;
        for (size_t j=0; j < n2; j++){
            traceno = this->sortmap[this->keymap[i].i0 + j];
            Fdata->seekg(Fdata->getStartofdata() + traceno*(n1+NHEAD2D)*sizeof(T));
            if(!this->getReciprocity()){
                Fdata->read(&scoords[j].x, 1);
                Fdata->read(&scoords[j].y, 1);
                Fdata->read(&gcoords[j].x, 1);
                Fdata->read(&gcoords[j].y, 1);
            }else{
                Fdata->read(&gcoords[j].x, 1);
                Fdata->read(&gcoords[j].y, 1);
                Fdata->read(&scoords[j].x, 1);
                Fdata->read(&scoords[j].y, 1);
            }
            Fdata->read(&data[j*n1], n1);
        }

        // Flag shot as running
        this->keymap[i].status = RUNNING;
        return gather;
    }else{
        // No shot available
        return nullptr;
    }
}

template<typename T>
std::shared_ptr<Data2D<T>> Sort<T>::get2DGather(size_t number)
{

    if(this->ngathers == 0 || this->ntraces == 0) rs_error("Sort::get2DGather: No sort map created.");
    if(number > ngathers-1) rs_error("Sort::get2DGather: Trying to get a gather with number that is larger than ngathers");

    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(this->datafile);
    if(status == FILE_ERR) rs_error("Sort::get2DGather: Error reading input data file: ", this->datafile);
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA2D) rs_error("Sort::get2DGather: Datafile must be of type Data2D.");
    //Get gather size information
    size_t n1 = Fdata->getN(1);
    T d1 = Fdata->getD(1);
    T o1 = Fdata->getO(1);
    size_t n2 = this->keymap[number].n;
    if(keymap[number].status == FINISHED)
    {
        rs_warning("Shot already processed.");
    }

    //Create gather
    std::shared_ptr<rockseis::Data2D<T>> gather;
    gather = std::make_shared<Data2D<T>>(n2,n1,d1,o1);
    Point2D<T> *scoords = (gather->getGeom())->getScoords();
    Point2D<T> *gcoords = (gather->getGeom())->getGcoords();
    T *data = gather->getData();
    size_t traceno;
    for (size_t j=0; j < n2; j++){
        traceno = this->sortmap[this->keymap[number].i0 + j];
        Fdata->seekg(Fdata->getStartofdata() + traceno*(n1+NHEAD2D)*sizeof(T));
        if(!this->getReciprocity()){
            Fdata->read(&scoords[j].x, 1);
            Fdata->read(&scoords[j].y, 1);
            Fdata->read(&gcoords[j].x, 1);
            Fdata->read(&gcoords[j].y, 1);
        }else{
            Fdata->read(&gcoords[j].x, 1);
            Fdata->read(&gcoords[j].y, 1);
            Fdata->read(&scoords[j].x, 1);
            Fdata->read(&scoords[j].y, 1);
        }
        Fdata->read(&data[j*n1], n1);
    }
    // Flag shot as running
    this->keymap[number].status = RUNNING;
    return gather;
}

template<typename T>
void Sort<T>::put2DGather(std::shared_ptr<Data2D<T>> data, size_t number)
{
    if(this->ngathers == 0 || this->ntraces == 0) rs_error("Sort::put2DGather: No sort map created.");
    if(number > ngathers-1) rs_error("Sort::put2DGather: Trying to put a gather with number that is larger than ngathers");

    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->append(data->getFile());
    if(status == FILE_ERR) rs_error("Sort::put2DGather: Error opening data file for appending: ", data->getFile());
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA2D) rs_error("Sort::put2DGather: Datafile must be of type Data2D.");
    //Get gather size information
    size_t n1 = Fdata->getN(1);
    T d1 = Fdata->getD(1);
    T o1 = Fdata->getO(1);
    size_t n2 = this->keymap[number].n;
    if(n1 != data->getNt()) rs_error("Sort::put2DGather: Number of samples in data and datafile mismatch.");
    if(d1 != data->getDt()) rs_error("Sort::put2DGather: Sampling interval in data and datafile mismatch.");
    if(o1 != data->getOt()) rs_error("Sort::put2DGather: Origin in data and datafile mismatch.");

    //Write gather
    Point2D<T> *scoords = (data->getGeom())->getScoords();
    Point2D<T> *gcoords = (data->getGeom())->getGcoords();
    T *tracedata = data->getData();
    size_t traceno;
    for (size_t j=0; j < n2; j++){
        traceno = this->sortmap[this->keymap[number].i0 + j];
        Fdata->seekp(Fdata->getStartofdata() + traceno*(n1+NHEAD2D)*sizeof(T));
        if(!this->getReciprocity()){
            Fdata->write(&scoords[j].x, 1);
            Fdata->write(&scoords[j].y, 1);
            Fdata->write(&gcoords[j].x, 1);
            Fdata->write(&gcoords[j].y, 1);
        }else{
            Fdata->write(&gcoords[j].x, 1);
            Fdata->write(&gcoords[j].y, 1);
            Fdata->write(&scoords[j].x, 1);
            Fdata->write(&scoords[j].y, 1);
        }
        Fdata->write(&tracedata[j*n1], n1);
    }

    if(Fdata->getFail()) rs_error("Sort::Put2DGather: Error writting gather to output file");
}

template<typename T>
void Sort<T>::put3DGather(std::shared_ptr<Data3D<T>> data, size_t number)
{
    if(this->ngathers == 0 || this->ntraces == 0) rs_error("Sort::put3DGather: No sort map created.");
    if(number > ngathers-1) rs_error("Sort::put3DGather: Trying to put a gather with number that is larger than ngathers");

    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->append(data->getFile());
    if(status == FILE_ERR) rs_error("Sort::put3DGather: Error reading input data file.");
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA3D) rs_error("Sort::put3DGather: Datafile must be of type Data3D.");
    //Get gather size information
    size_t n1 = Fdata->getN(1);
    T d1 = Fdata->getD(1);
    T o1 = Fdata->getO(1);
    size_t n2 = this->keymap[number].n;
    if(n1 != data->getNt()) rs_error("Sort::put3DGather: Number of samples in data and datafile mismatch.");
    if(d1 != data->getDt()) rs_error("Sort::put3DGather: Sampling interval in data and datafile mismatch.");
    if(o1 != data->getOt()) rs_error("Sort::put3DGather: Origin in data and datafile mismatch.");

    //Write gather
    Point3D<T> *scoords = (data->getGeom())->getScoords();
    Point3D<T> *gcoords = (data->getGeom())->getGcoords();
    T *tracedata = data->getData();
    size_t traceno;
    for (size_t j=0; j < n2; j++){
        traceno = this->sortmap[this->keymap[number].i0 + j];
        Fdata->seekp(Fdata->getStartofdata() + traceno*(n1+NHEAD3D)*sizeof(T));
        if(!this->getReciprocity()){
            Fdata->write(&scoords[j].x, 1);
            Fdata->write(&scoords[j].y, 1);
            Fdata->write(&scoords[j].z, 1);
            Fdata->write(&gcoords[j].x, 1);
            Fdata->write(&gcoords[j].y, 1);
            Fdata->write(&gcoords[j].z, 1);
        }else{
            Fdata->write(&gcoords[j].x, 1);
            Fdata->write(&gcoords[j].y, 1);
            Fdata->write(&gcoords[j].z, 1);
            Fdata->write(&scoords[j].x, 1);
            Fdata->write(&scoords[j].y, 1);
            Fdata->write(&scoords[j].z, 1);
        }
        Fdata->write(&tracedata[j*n1], n1);
    }

    if(Fdata->getFail()) rs_error("Sort::Put3DGather: Error writting gather to output file");
}

template<typename T>
std::shared_ptr<Data3D<T>> Sort<T>::get3DGather()
{
    if(this->ngathers == 0 || this->ntraces == 0) rs_error("No sort map created.");

    int i;
    for(i=0; i < this->ngathers; i++)
    {
        if(keymap[i].status == NOT_STARTED || keymap[i].status == FAILED)
        {
            break;
        }
    }
    if(i < this->ngathers) 
    {
        // Found a shot
        bool status;
        std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
        if(strcmp(this->datafile.c_str(), "stdin")){
            status = Fdata->input(this->datafile);
        }else{
            status = Fdata->input();
        }
        if(status == FILE_ERR) rs_error("Sort::get3DGather: Error reading input data file: ", this->datafile);
        rs_datatype datatype = Fdata->getType(); 
        if(datatype != DATA3D) rs_error("Sort::get3DGather: Datafile must be of type Data3D.");
        //Get gather size information
        size_t n1 = Fdata->getN(1);
        T d1 = Fdata->getD(1);
        T o1 = Fdata->getO(1);
        size_t n2 = this->keymap[i].n;

        //Create gather
	    std::shared_ptr<rockseis::Data3D<T>> gather;
        gather = std::make_shared<Data3D<T>>(n2,n1,d1,o1);
        Point3D<T> *scoords = (gather->getGeom())->getScoords();
        Point3D<T> *gcoords = (gather->getGeom())->getGcoords();
        T *data = gather->getData();
        size_t traceno;
        for (size_t j=0; j < n2; j++){
            traceno = this->sortmap[this->keymap[i].i0 + j];
            Fdata->seekg(Fdata->getStartofdata() + traceno*(n1+NHEAD3D)*sizeof(T));
            if(!this->getReciprocity()){
                Fdata->read(&scoords[j].x, 1);
                Fdata->read(&scoords[j].y, 1);
                Fdata->read(&scoords[j].z, 1);
                Fdata->read(&gcoords[j].x, 1);
                Fdata->read(&gcoords[j].y, 1);
                Fdata->read(&gcoords[j].z, 1);
            }else{
                Fdata->read(&gcoords[j].x, 1);
                Fdata->read(&gcoords[j].y, 1);
                Fdata->read(&gcoords[j].z, 1);
                Fdata->read(&scoords[j].x, 1);
                Fdata->read(&scoords[j].y, 1);
                Fdata->read(&scoords[j].z, 1);
            }
            Fdata->read(&data[j*n1], n1);
        }
        // Flag shot as running
        this->keymap[i].status = RUNNING;
        return gather;
    }else{
        // No shot available
        return nullptr;
    }
}

template<typename T>
std::shared_ptr<Data3D<T>> Sort<T>::get3DGather(size_t number)
{
    if(this->ngathers == 0 || this->ntraces == 0) rs_error("No sort map created.");

    if(number > ngathers-1) rs_error("Sort::get3DGather: Trying to get a gather with number that is larger than ngathers");


    // Found a shot
    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(this->datafile);
    if(status == FILE_ERR) rs_error("Sort::get3DGather: Error reading input data file: ", this->datafile);
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA3D) rs_error("Sort::get3DGather: Datafile must be of type Data3D.");
    //Get gather size information
    size_t n1 = Fdata->getN(1);
    T d1 = Fdata->getD(1);
    T o1 = Fdata->getO(1);
    size_t n2 = this->keymap[number].n;

    //Create gather
    std::shared_ptr<rockseis::Data3D<T>> gather;
    gather = std::make_shared<Data3D<T>>(n2,n1,d1,o1);
    Point3D<T> *scoords = (gather->getGeom())->getScoords();
    Point3D<T> *gcoords = (gather->getGeom())->getGcoords();
    T *data = gather->getData();
    size_t traceno;
    for (size_t j=0; j < n2; j++){
        traceno = this->sortmap[this->keymap[number].i0 + j];
        Fdata->seekg(Fdata->getStartofdata() + traceno*(n1+NHEAD3D)*sizeof(T));
        if(!this->getReciprocity()){
            Fdata->read(&scoords[j].x, 1);
            Fdata->read(&scoords[j].y, 1);
            Fdata->read(&scoords[j].z, 1);
            Fdata->read(&gcoords[j].x, 1);
            Fdata->read(&gcoords[j].y, 1);
            Fdata->read(&gcoords[j].z, 1);
        }else{
            Fdata->read(&gcoords[j].x, 1);
            Fdata->read(&gcoords[j].y, 1);
            Fdata->read(&gcoords[j].z, 1);
            Fdata->read(&scoords[j].x, 1);
            Fdata->read(&scoords[j].y, 1);
            Fdata->read(&scoords[j].z, 1);
        }
        Fdata->read(&data[j*n1], n1);
    }
    // Flag shot as running
    this->keymap[number].status = RUNNING;
    return gather;
}


template<typename T>
void Sort<T>::readKeymap(){
    if(this->kmapfile.empty()) rs_error("Sort::readKeymap: kmapfile not set.");
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    if(Fin->input(this->kmapfile) == FILE_ERR) rs_error("Sort::readKeymap: Error opening file: ", this->kmapfile);
    if(Fin->getType() != KEYMAP) rs_error("Sort::readKeymap: Key map file is not of correct type (KEYMAP)");
    this->ngathers = Fin->getN(1) - 1;
    free(keymap);
    keymap = (key *) calloc(this->ngathers, sizeof(key));
    int status;
    size_t reciprocityflag;
    for (size_t i=0; i < this->ngathers; i++)
    {
        Fin->read(&(this->keymap[i].i0), 1);
        Fin->read(&(this->keymap[i].n), 1);
        Fin->read(&status, 1);
        this->keymap[i].status = static_cast<rockseis::rs_status>(status);
    }
    // Read one more value with sort key and reciprocity flag
    Fin->read(&(reciprocityflag), 1);
    Fin->read(&(reciprocityflag), 1);
    if(reciprocityflag)
        reciprocity = true;
    else{
        reciprocity = false;
    }
    Fin->read(&status, 1);
    sortkey = static_cast<rs_key>(status);
    Fin->close();
    sortset = true;
}

template<typename T>
size_t Sort<T>::getMaxtraces(){
    if(!this->sortset) rs_error("Sort::getMaxtraces: sort map is not created.");

    size_t max = 0;
    for (size_t i=0; i < this->ngathers; i++)
    {
        if(this->keymap[i].n > max) max = this->keymap[i].n;
    }
    return max;
}

template<typename T>
size_t Sort<T>::getMintraces(){
    if(!this->sortset) rs_error("Sort::getMintraces: sort map is not created.");

    size_t min = 0;
    for (size_t i=0; i < this->ngathers; i++)
    {
        if(this->keymap[i].n < min) min = this->keymap[i].n;
    }
    return min;
}

template<typename T>
void Sort<T>::writeKeymap(){
    if(ngathers == 0) rs_error("Sort::writeKeymap: Key map was not created, ngathers = 0");
    size_t n1 = this->ngathers;

    if(this->kmapfile.empty()) rs_error("Sort::writeKeymap: kmapfile not set.");

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    Fout->output(this->kmapfile);
    Fout->setN(1,n1+1);
    Fout->setD(1,1);
    Fout->setData_format(2*sizeof(size_t) + sizeof(int));
    Fout->setType(KEYMAP);
    Fout->writeHeader();
    Fout->seekp(Fout->getStartofdata());
    int status;
    size_t reciprocityflag;
    for (size_t i=0; i < n1; i++)
    {
        Fout->write(&(this->keymap[i].i0), 1);
        Fout->write(&(this->keymap[i].n), 1);
        status = static_cast<int>(this->keymap[i].status);
        Fout->write(&status, 1);
    }
    // Write one more value with sort key and reciprocity flag
    if(this->getReciprocity())
        reciprocityflag = 1;
    else{
        reciprocityflag = 0;
    }
    Fout->write(&(reciprocityflag), 1);
    Fout->write(&(reciprocityflag), 1);
    status = static_cast<int>(this->getSortkey());
    Fout->write(&status, 1);
    Fout->close();
}

template<typename T>
void Sort<T>::readSortmap(){
    if(this->smapfile.empty()) rs_error("Sort::readSortmap: smapfile not set.");
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    if(Fin->input(this->smapfile) == FILE_ERR) rs_error("Sort::readSortmap: Error opening file: ", this->smapfile);
    if(Fin->getType() != SORTMAP) rs_error("Sort::readSortmap: Sort map file is not of correct type (SORTMAP)");
    this->ntraces = Fin->getN(1);
    free(sortmap);
    sortmap = (size_t *) calloc(this->ntraces, sizeof(size_t));
    for (size_t i=0; i < this->ntraces; i++)
    {
        Fin->read(&(this->sortmap[i]), 1);
    }
}

template<typename T>
void Sort<T>::writeSortmap(){
    if(ngathers == 0) rs_error("Sort::writeSortmap: Sort map was not created, ngathers = 0");
    size_t n1 = this->ntraces;

    if(this->smapfile.empty()) rs_error("Sort::writeSortmap: smapfile not set.");

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    Fout->output(this->smapfile);
    Fout->setN(1,n1);
    Fout->setD(1,1);
    Fout->setData_format(sizeof(size_t));
    Fout->setType(SORTMAP);
    Fout->writeHeader();
    Fout->seekp(Fout->getStartofdata());
    for (size_t i=0; i < n1; i++)
    {
        Fout->write(&(this->sortmap[i]), 1);
    }
}


template<typename T>
void Sort<T>::createEmptydataset(std::string filename, size_t n1, T d1, T o1){
    if(this->ngathers == 0 || this->ntraces == 0) rs_error("No sort map created.");

    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(this->datafile);
    if(status == FILE_ERR) rs_error("Sort::createEmptydata: Error reading input survey data file: ", this->datafile);
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA2D && datatype != DATA3D ) rs_error("Sort::createEmptydata: Datafile must be of type Data2D or Data3D.");

    //Get number of traces
    size_t n2 = this->ntraces;
    size_t nt = Fdata->getN(1);
    std::shared_ptr<rockseis::Data2D<T>> data2d;
    std::shared_ptr<rockseis::Data3D<T>> data3d;
    if(datatype == DATA2D){
        data2d = std::make_shared<Data2D<T>>(1,n1,d1,o1);
        data2d->setFile(filename);
        status = data2d->open("o");
        if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

        Point2D<T> *scoords = (data2d->getGeom())->getScoords();
        Point2D<T> *gcoords = (data2d->getGeom())->getGcoords();
        for (size_t j=0; j < n2; j++){
            Fdata->seekg(Fdata->getStartofdata() + j*(nt+NHEAD2D)*sizeof(T));
            Fdata->read(&scoords[0].x, 1);
            Fdata->read(&scoords[0].y, 1);
            Fdata->read(&gcoords[0].x, 1);
            Fdata->read(&gcoords[0].y, 1);
            data2d->writeTraces();
        }
        data2d->close();
        if(Fdata->getFail()) rs_error("Sort::createEmptydataset: Error reading from survey file");
    }else{
        data3d = std::make_shared<Data3D<T>>(1,n1,d1,o1);
        data3d->setFile(filename);
        status = data3d->open("o");
        if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

        Point3D<T> *scoords = (data3d->getGeom())->getScoords();
        Point3D<T> *gcoords = (data3d->getGeom())->getGcoords();
        for (size_t j=0; j < n2; j++){
            Fdata->seekg(Fdata->getStartofdata() + j*(nt+NHEAD3D)*sizeof(T));
            Fdata->read(&scoords[0].x, 1);
            Fdata->read(&scoords[0].y, 1);
            Fdata->read(&scoords[0].z, 1);
            Fdata->read(&gcoords[0].x, 1);
            Fdata->read(&gcoords[0].y, 1);
            Fdata->read(&gcoords[0].z, 1);
            data3d->writeTraces();
        }
        data3d->close();
        if(Fdata->getFail()) rs_error("Sort::createEmptydataset: Error reading from survey file");
    }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Sort<float>;
template class Sort<double>;

}
