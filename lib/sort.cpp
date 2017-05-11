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

	if( pa->x < pb->x) return -1;
	if( pa->x > pb->x) return 1;

	if( pa->y < pb->y) return -1;
	if( pa->y > pb->y) return 1;

	if( pa->z < pb->z) return -1;
	if( pa->z > pb->z) return 1;

	return 0;
}

//Sort with distance from origin as primary key and y as secondary key
int sort_sr_1(void const *a, void const *b)
{
	position_t *pa, *pb;

	pa = (position_t *) a;
	pb = (position_t *) b;

	if((pa->x*pa->x) + (pa->y*pa->y) < (pb->x*pb->x) + (pb->y*pb->y)) return -1;
	if((pa->x*pa->x) + (pa->y*pa->y) > (pb->x*pb->x) + (pb->y*pb->y)) return 1;

	if( pa->y < pb->y) return -1;
	if( pa->y > pb->y) return 1;

	if( pa->z < pb->z) return -1;
	if( pa->z > pb->z) return 1;

	return 0;
}

template<typename T>
bool Sort<T>::createSort(std::string filename, rs_key _sortkey, T dx, T dy)
{
    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(filename);
    if(status == FILE_ERR) rs_error("Sort::createSort: Error reading input file.");
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
                            positions[i2].y = fbuffer[1];
                            positions[i2].z = fbuffer[2];
                        }else{
                            positions[i2].x = fbuffer[0];
                            positions[i2].y = 0.0;
                            positions[i2].z = fbuffer[1];
                        }
                        break;
                    case RECEIVER:
                        if(datatype == DATA3D) {
                            positions[i2].x = fbuffer[3];
                            positions[i2].y = fbuffer[4];
                            positions[i2].z = fbuffer[5];
                        }else{
                            positions[i2].x = fbuffer[2];
                            positions[i2].y = 0.0;
                            positions[i2].z = fbuffer[3];
                        }
                        break;
                    case CMP:
                        if(datatype == DATA3D) {
                            positions[i2].x = 0.5*(fbuffer[0]+fbuffer[3]);
                            positions[i2].y = 0.5*(fbuffer[1]+fbuffer[4]);
                        }else{
                            positions[i2].x = 0.5*(fbuffer[0]+fbuffer[2]);
                            positions[i2].y = 0.0;
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
                            positions[i2].y = dbuffer[1];
                            positions[i2].z = dbuffer[2];
                        }else{
                            positions[i2].x = dbuffer[0];
                            positions[i2].y = 0.0;
                            positions[i2].z = dbuffer[1];
                        }
                        break;
                    case RECEIVER:
                        if(datatype == DATA3D) {
                            positions[i2].x = dbuffer[3];
                            positions[i2].y = dbuffer[4];
                            positions[i2].z = dbuffer[5];
                        }else{
                            positions[i2].x = dbuffer[2];
                            positions[i2].y = 0.0;
                            positions[i2].z = dbuffer[3];
                        }
                        break;
                    case CMP:
                        if(datatype == DATA3D) {
                            positions[i2].x = 0.5*(dbuffer[0]+dbuffer[3]);
                            positions[i2].y = 0.5*(dbuffer[1]+dbuffer[4]);
                        }else{
                            positions[i2].x = 0.5*(dbuffer[0]+dbuffer[2]);
                            positions[i2].y = 0.0;
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
        status = Fdata->input(this->datafile);
        if(status == FILE_ERR) rs_error("Sort::get2DGather: Error reading input data file.");
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
            Fdata->read(&scoords[j].x, 1);
            Fdata->read(&scoords[j].y, 1);
            Fdata->read(&gcoords[j].x, 1);
            Fdata->read(&gcoords[j].y, 1);
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
    if(status == FILE_ERR) rs_error("Sort::get2DGather: Error reading input data file.");
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
        Fdata->read(&scoords[j].x, 1);
        Fdata->read(&scoords[j].y, 1);
        Fdata->read(&gcoords[j].x, 1);
        Fdata->read(&gcoords[j].y, 1);
        Fdata->read(&data[j*n1], n1);
    }
    // Flag shot as running
    this->keymap[number].status = RUNNING;
    return gather;
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
        status = Fdata->input(this->datafile);
        if(status == FILE_ERR) rs_error("Sort::get3DGather: Error reading input data file.");
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
            Fdata->read(&scoords[j].x, 1);
            Fdata->read(&scoords[j].y, 1);
            Fdata->read(&scoords[j].z, 1);
            Fdata->read(&gcoords[j].x, 1);
            Fdata->read(&gcoords[j].y, 1);
            Fdata->read(&gcoords[j].z, 1);
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

    if(number > ngathers-1) rs_error("Sort::get2DGather: Trying to get a gather with number that is larger than ngathers");


    // Found a shot
    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(this->datafile);
    if(status == FILE_ERR) rs_error("Sort::get3DGather: Error reading input data file.");
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
        Fdata->read(&scoords[j].x, 1);
        Fdata->read(&scoords[j].y, 1);
        Fdata->read(&scoords[j].z, 1);
        Fdata->read(&gcoords[j].x, 1);
        Fdata->read(&gcoords[j].y, 1);
        Fdata->read(&gcoords[j].z, 1);
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
    this->ngathers = Fin->getN(1);
    free(keymap);
    keymap = (key *) calloc(this->ngathers, sizeof(key));
    int status;
    for (size_t i=0; i < this->ngathers; i++)
    {
        Fin->read(&(this->keymap[i].i0), 1);
        Fin->read(&(this->keymap[i].n), 1);
        Fin->read(&status, 1);
        this->keymap[i].status = static_cast<rockseis::rs_status>(status);
    }
}

template<typename T>
void Sort<T>::writeKeymap(){
    if(ngathers == 0) rs_error("Sort::writeKeymap: Key map was not created, ngathers = 0");
    size_t n1 = this->ngathers;

    if(this->kmapfile.empty()) rs_error("Sort::writeKeymap: kmapfile not set.");

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    Fout->output(this->kmapfile);
    Fout->setN(1,n1);
    Fout->setD(1,1);
    Fout->setData_format(2*sizeof(size_t) + sizeof(int));
    Fout->setType(KEYMAP);
    Fout->writeHeader();
    Fout->seekp(Fout->getStartofdata());
    int status;
    for (size_t i=0; i < n1; i++)
    {
        Fout->write(&(this->keymap[i].i0), 1);
        Fout->write(&(this->keymap[i].n), 1);
        status = static_cast<int>(this->keymap[i].status);
        Fout->write(&status, 1);
    }
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
    if(status == FILE_ERR) rs_error("Sort::createEmptydata: Error reading input survey data file.");
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA2D && datatype != DATA3D ) rs_error("Sort::createEmptydata: Datafile must be of type Data2D or Data3D.");

    //Get number of traces
    size_t n2 = this->ntraces;
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
            Fdata->seekg(Fdata->getStartofdata() + j*(n1+NHEAD2D)*sizeof(T));
            Fdata->read(&scoords[0].x, 1);
            Fdata->read(&scoords[0].y, 1);
            Fdata->read(&gcoords[0].x, 1);
            Fdata->read(&gcoords[0].y, 1);
            data2d->writeTraces();
        }
        data2d->close();
    }else{
        data3d = std::make_shared<Data3D<T>>(1,n1,d1,o1);
        data3d->setFile(filename);
        status = data3d->open("o");
        if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

        Point3D<T> *scoords = (data3d->getGeom())->getScoords();
        Point3D<T> *gcoords = (data3d->getGeom())->getGcoords();
        for (size_t j=0; j < n2; j++){
            Fdata->seekg(Fdata->getStartofdata() + j*(n1+NHEAD3D)*sizeof(T));
            Fdata->read(&scoords[0].x, 1);
            Fdata->read(&scoords[0].y, 1);
            Fdata->read(&scoords[0].z, 1);
            Fdata->read(&gcoords[0].x, 1);
            Fdata->read(&gcoords[0].y, 1);
            Fdata->read(&gcoords[0].z, 1);
            data3d->writeTraces();
        }
        data3d->close();
    }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Sort<float>;
template class Sort<double>;

}
