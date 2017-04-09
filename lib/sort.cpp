#include "sort.h"

namespace rockseis {
// constructor
template<typename T>
Sort<T>::Sort()
{
    nensembles = 0;
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
int sort_xy_0(void const *a, void const *b)
{
	position_t *pa, *pb;

	pa = (position_t *) a;
	pb = (position_t *) b;

	if( pa->x < pb->x) return -1;
	if( pa->x > pb->x) return 1;

	if( pa->y < pb->y) return -1;
	if( pa->y > pb->y) return 1;

	return 0;
}

//Sort with distance from origin as primary key and y as secondary key
int sort_xy_1(void const *a, void const *b)
{
	position_t *pa, *pb;

	pa = (position_t *) a;
	pb = (position_t *) b;

	if((pa->x*pa->x) + (pa->y*pa->y) < (pb->x*pb->x) + (pb->y*pb->y)) return -1;
	if((pa->x*pa->x) + (pa->y*pa->y) > (pb->x*pb->x) + (pb->y*pb->y)) return 1;

	if( pa->y < pb->y) return -1;
	if( pa->y > pb->y) return 1;

	return 0;
}

template<typename T>
bool Sort<T>::createSort(std::string filename, rs_key _sortkey)
{
    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(filename);
    if(status == FILE_ERR) rs_error("Sort::createSort: Error reading input file.");
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA2D && datatype != DATA3D) rs_error("Sort::createSort: Only DATA2D and DATA3D types are supported.");

    sortkey = _sortkey; 
    if(sortkey != SOURCE && sortkey != RECEIVER) rs_error("Sort::createSort: Only SOURCE and RECEIVER sort keys are supported.");
    size_t n1 = Fdata->getN(1);
    size_t n2 = Fdata->getN(2);

    sortmap = (size_t *) calloc(n2, sizeof(size_t));

    int nh = Fdata->getNheader();
    int dp = Fdata->getData_format();
    int hp = Fdata->getHeader_format();
    float fbuffer[6];
    double dbuffer[6];
    position_t *positions, *ensembles;

	positions = (position_t *) calloc(n2,sizeof(position_t));
	ensembles = (position_t *) calloc(n2,sizeof(position_t));

    nensembles = 0;
    for (int i2=0; i2 < n2; i2++){
        //Read coordinates
        Fdata->seekg(i2*(n1*dp + nh*hp));
        switch (hp){
            case sizeof(float):
                Fdata->read(fbuffer, nh);
                switch (sortkey){
                    case SOURCE:
                        positions[i2].x = fbuffer[0];
                        positions[i2].y = fbuffer[1];
                        break;
                    case RECEIVER:
                        positions[i2].x = fbuffer[2];
                        positions[i2].y = fbuffer[3];
                        break;
                    default:
                        break;
                }
                break;
            case sizeof(double):
                Fdata->read(dbuffer, nh);
                switch (sortkey){
                    case SOURCE:
                        positions[i2].x = dbuffer[0];
                        positions[i2].y = dbuffer[1];
                        break;
                    case RECEIVER:
                        positions[i2].x = dbuffer[2];
                        positions[i2].y = dbuffer[3];
                        break;
                    default:
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
	qsort(positions, n2, sizeof(position_t), sort_xy_0);

    nensembles = 1;
    int ntraces = 1;
    double oldx, oldy;
    size_t oldind;

	/* Initialize the ensemble list and map */
	oldx = positions[0].x;
	oldy = positions[0].y;
	ensembles[0].x=oldx;
	ensembles[0].y=oldy;
	oldind = 0;

	// Making a map of the traces 
	size_t *map = (size_t *) calloc(2*n2,sizeof(size_t));

	for (int i2=1; i2 < n2; i2++) {
		// Get coordinates
		dbuffer[0] = positions[i2].x;
		dbuffer[1] = positions[i2].y;
		// If ensemble position is not in the list then add it to the list
		if((dbuffer[0] == oldx) && (dbuffer[1] == oldy)){
			ntraces++;
		}else{
			ensembles[nensembles].x=dbuffer[0];
			ensembles[nensembles].y=dbuffer[1];
			map[2*(nensembles-1)]=oldind;
			map[2*(nensembles-1)+1]=ntraces;
			nensembles++;
			ntraces=1;
			oldx=dbuffer[0];
			oldy=dbuffer[1];
			oldind=i2;
		}
    }
    // Taking the last ensemble into consideration
	map[2*(nensembles-1)]=oldind;
	map[2*(nensembles-1)+1]=ntraces;

    // Create keymap and sortmap
    free(keymap); free(sortmap);
    keymap = (key *) calloc(nensembles, sizeof(key));
    sortmap = (size_t *) calloc(n2, sizeof(size_t));

	for (int i2=0; i2 < nensembles; i2++) {
        keymap[i2].i0 = map[2*i2];
        keymap[i2].n = map[2*i2 + 1];
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

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Sort<float>;
template class Sort<double>;

}
