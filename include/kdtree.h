#ifndef _kdtree_h
#define _kdtree_h


struct kdtree {
	int dim;
	struct kdnode *root;
	struct kdhyperrect *rect;
	void (*destr)(void*);
};


struct kdres {
	struct kdtree *tree;
	struct res_node *rlist, *riter;
	int size;
};


struct kdtree *kd_create(int k);
/*<create a kd-tree for "k"-dimensional data>*/


void kd_clear(struct kdtree *tree);
/*<remove all the elements from the tree >*/


void kd_free(struct kdtree *tree);
/*<free the struct kdtree>*/


void kd_data_destructor(struct kdtree *tree, void (*destr)(void*));
/*<if called with non-null 2nd argument, the function provided will be called on data pointers (see kd_insert) when nodes are to be removed from the tree.>*/


int kd_insert(struct kdtree *tree, const double *pos, int data);
/*<insert a node, specifying its position, and optional data >*/


int kd_insertf(struct kdtree *tree, const float *pos, int data);
/*<insert a node, specifying its position, and optional data >*/


int kd_insert3(struct kdtree *tree, double x, double y, double z, int data);
/*<insert a node, specifying its position, and optional data >*/


int kd_insert3f(struct kdtree *tree, float x, float y, float z, int data);
/*<insert a node, specifying its position, and optional data >*/


void kd_res_free(struct kdres *rset);
/*<frees a result set returned by kd_nearest_range() >*/


void kd_res_rewind(struct kdres *rset);
/*<rewinds the result set iterator >*/


struct kdres *kd_nearest_range(struct kdtree *kd, const double *pos, double range);
/*<Find any nearest nodes from a given point within a range.>*/


struct kdres *kd_nearest_rangef(struct kdtree *kd, const float *pos, float range);
/*<Find any nearest nodes from a given point within a range.>*/


struct kdres *kd_nearest_range3(struct kdtree *tree, double x, double y, double z, double range);
/*<Find any nearest nodes from a given point within a range.>*/


struct kdres *kd_nearest_range3f(struct kdtree *tree, float x, float y, float z, float range);
/*<Find any nearest nodes from a given point within a range.>*/


int kd_res_size(struct kdres *set);
/*<returns the size of the result set (in elements) >*/


int kd_res_end(struct kdres *rset);
/*<returns non-zero if the set iterator reached the end after the last element>*/


int kd_res_next(struct kdres *rset);
/*<advances the result set iterator, returns non-zero on success, zero if there are no more elements in the result set.>*/


int kd_res_item(struct kdres *rset, double *pos);
/*<returns the data pointer (can be null) of the current result set item and optionally sets its position to the pointers(s) if not null.>*/


int kd_res_itemf(struct kdres *rset, float *pos);
/*<returns the data pointer (can be null) of the current result set item and optionally sets its position to the pointers(s) if not null.>*/


int kd_res_item3(struct kdres *rset, double *x, double *y, double *z);
/*<returns the data pointer (can be null) of the current result set item and optionally sets its position to the pointers(s) if not null.>*/


int kd_res_item3f(struct kdres *rset, float *x, float *y, float *z);
/*<returns the data pointer (can be null) of the current result set item and optionally sets its position to the pointers(s) if not null.>*/


int kd_res_item_data(struct kdres *set);
/*<equivalent to kd_res_item(set, 0) >*/


int kd_find_nearest4(struct kdtree *tree, double *pos, double R, double *min_pos, int *min_index);
/*<Find nearest trace within the range specified by R>*/


int kd_find_nearest4f(struct kdtree *tree, float *pos, float R, float *min_pos, int *min_index);
/*<Find nearest trace within the range specified by R>*/

int kd_find_nearest2f(struct kdtree *tree, float *pos, float R, float *min_pos, int *min_index);
/*<Find nearest trace within the range specified by R>*/

#endif
