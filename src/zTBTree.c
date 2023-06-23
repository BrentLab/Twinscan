#include "zTBTree.h"

/**********************************************************************\
  Traceback Tree

  This is a tree structure which can store all live traceback paths 
  from a trellis more efficiently than retaining the full trellis.

\**********************************************************************/

coor_t any_new_cpoint = 0;
static int nodeid = 0;
zTBTreeNode* watchme; 
zTBTreeNode* watchme2;

static void* zCreateTBTreeNode(){
	zTBTreeNode* t = zMalloc(sizeof(zTBTreeNode),"zCreateTBTreeNode t");
	t->pos = 0;
	t->state = -1;
	t->score = MIN_SCORE;
	t->frame_data = UNDEFINED_FRAME;
	t->phase = Phase0;
	t->strand = UNDEFINED_STRAND;
	t->parent = NULL;
	t->child = NULL;
	t->lsib = NULL;
	t->rsib = NULL;
	t->children = 0;
	t->id = nodeid++;
	t->lock = 0;
	t->alive = 1;
	/*EVAN/
	for(i = 0;i <= 10; i++){
		t->ng[i] =0;
		}*/
	return (void*)t;
}

static void zFreeTBTreeNode(void* v){
	zTBTreeNode* n = (zTBTreeNode*)v;
	zFree(n);
}

static void zResetTBTreeNode(void* v){
	zTBTreeNode* n = (zTBTreeNode*)v;
	n->pos = 0;
	n->state = -1;
	n->score = MIN_SCORE;
	n->frame_data = UNDEFINED_FRAME;
	n->phase = Phase0;
	n->strand = UNDEFINED_STRAND;
	n->parent = NULL;
	n->child = NULL;
	n->lsib = NULL;
	n->rsib = NULL;
	n->children = 0;
	n->lock = 0;
	n->alive = 1;
	/*EVAN*
	for(i = 0;i <= 10; i++){
		n->ng[i] =0;
		}*/
}

zTBTreeNode* zGetTBTreeNode(zTBTree *t){
	zTBTreeNode* n;
	t->size++;
	if(t->dead_nodes == NULL){
		/* if there are no dead node allocate a new node */
		n = zCreateTBTreeNode();
	}
	else{
		/* otherwise return the dead node at the start of the list */
		n = t->dead_nodes;
		t->dead_nodes = t->dead_nodes->child;
		if(n->alive != 0){
			zDie("");
		}
		zResetTBTreeNode(n);
	}
	n->tree = t;
	return n;
}

static int TBTREE_ID = 0;

void zInitTBTree(zTBTree *t){
	t->dead_nodes = NULL;
	t->size = 0;
	t->root = zGetTBTreeNode(t);
	t->cpoint = t->root;
	t->new_cpoint = 0;
	t->id = TBTREE_ID++;
	t->check_count = 1;
}

static void zFreeTBTreeHelp(zTBTreeNode* n){
	zTBTreeNode* c = n->child;
	zTBTreeNode* c2;
	while(c != NULL){
		c2 = c;
		c = c->rsib;
		zFreeTBTreeHelp(c2);
	}
	zFreeTBTreeNode(n);
}

void zFreeTBTree(zTBTree *t){
	zTBTreeNode* n;
	zFreeTBTreeHelp(t->root);
	while(t->dead_nodes != NULL){
		n = t->dead_nodes;
		t->dead_nodes = t->dead_nodes->child;
		zFreeTBTreeNode(n);
	}
}

static void zReleaseTBTreeNode(zTBTree *t,zTBTreeNode* n){
	/* clear n and prepend it to the dead_nodes list */
	zResetTBTreeNode(n);
	n->alive = 0;
	n->child = t->dead_nodes;
	t->dead_nodes = n;
	t->size--;
}

static void zResetTBTreeHelp(zTBTree*,zTBTreeNode*);

void zClearTBTreeNodeChildren(zTBTree* t, zTBTreeNode* n){
	zTBTreeNode* c;
	if(n->children > 0){
		c = n->child;
		n->children = 0;
		n->child = NULL;
		while(c != NULL){
			n = c;
			c = c->rsib;
			zResetTBTreeHelp(t,n);
		}
	}
}

static void zResetTBTreeHelp(zTBTree* t, zTBTreeNode* n){
	zClearTBTreeNodeChildren(t,n);
	zReleaseTBTreeNode(t,n);
}

void zResetTBTree(zTBTree* tree){
	zResetTBTreeHelp(tree,tree->root);
	tree->root = zGetTBTreeNode(tree);
	tree->cpoint = tree->root;
	tree->new_cpoint = 0;
}

bool zTBTreeCheckNewCpoint(zTBTree* tree){
	return !(tree->new_cpoint == 0);
}

bool zTBTreeClearNewCpoint(zTBTree* tree){
	if(tree->new_cpoint){
		tree->new_cpoint = 0;
		return true;
	}
	else{
		return false;
	}
}

/*EVAN
static void zTBTreeVerifyStructureHelp(zTBTree* t, zTBTreeNode* n){
	zTBTreeNode* c = n->child;
	int count = 0;
	if(n->tree != t){
		zDie("");
	}
	if(n->alive == 0){
		zDie("");
	}
	if(n->children == 0){
		if(c != NULL){
			zDie("");
		}
	}
	else if(n->children < 0){
		zDie("");
	}
	else{
		if(c == NULL){
			zDie("");
		}
		if(c->lsib != NULL){
			zDie("");
		}
		while(c != NULL){
			count++;
			if(c->parent != n){
				zDie("");
			}
			if(c->rsib != NULL){
				if(c->rsib->lsib != c){
					zDie("");
				}
			}
			zTBTreeVerifyStructureHelp(t,c);
			c = c->rsib;
		}
		if(count != n->children){
			zDie("");
		}
	}
}

static void zTBTreeVerifyStructure(zTBTree* t){
	zTBTreeNode* n = t->dead_nodes;
	t->check_count++;
	while(n != NULL){
		if(n->tree != t){
			zDie("");
		}
		if(n->alive != 0){
			zDie("");
		}
		n = n->child;
	}
	zTBTreeVerifyStructureHelp(t,t->root);
}

static bool zTBTreeNodeCheckChildren(zTBTreeNode* n){
	zTBTreeNode* c;
	int count = 0;
	bool val = true;
	if(n->pos == 0 && n->state == -1){
		return true;
	}
	if(n->child == NULL){
		if(n->children != 0){
			val = false;
		}
	}
	else{
		count = 1;
		c = n->child;
		if(c->parent != n){
			zDie("WTF\n");
		}
		while(c->rsib != NULL){
			c = c->rsib;
			count++;
			if(c->parent != n){
				zDie("WTF\n");
			}
		}
		if(count != n->children){
			val = false;
		}
	}
	if(val == false){
		zDie("FINDME\n");
	}
	return val;
}
*/

/* this assumes parent and 1 child. child must be same state as 
   this node or we are throwing away info */
void zReleaseRedundantTBTreeNode(zTBTree* t,zTBTreeNode* n){
	zTBTreeNode* p = n->parent;
	zTBTreeNode* c = n->child;
	if(n->lock == 0){
		c->parent = p;
		/*EVAN heuristic
		for(i = 0; i <= 10; i++){
			if(n->ng[i]){
				c->ng[i] = 1;
			}
			}*/
		if(p->child == n){
			p->child = c;
		}
		c->lsib = n->lsib;
		c->rsib = n->rsib;
		if(c->rsib != NULL){
			c->rsib->lsib = c;
		}
		if(c->lsib != NULL){
			c->lsib->rsib = c;
		}
	}
	if(n == t->cpoint){
		t->cpoint = n->child;
		while(t->cpoint->children == 1){
			t->cpoint = t->cpoint->child;
		}
		t->new_cpoint = 1;
		any_new_cpoint = t->cpoint->pos;
	}
	if(n->lock == 0){
		zReleaseTBTreeNode(t,n);
	}
}

/* this assumes that the node has no children and a parent */
void zReleaseDeadTBTreeNode(zTBTree* t, zTBTreeNode* n){
	zTBTreeNode* dead;
	/* walk tree toward root removing new dead nodes */ 
	while(n->parent->children == 1){
		n = n->parent;
		zReleaseTBTreeNode(t,n->child);
		n->children = 0;
		n->child = NULL;
	}
	/* remove current node */
	dead = n;
	n = n->parent;
	if(n->children <= 1){
		zDie("hwhap\n");
	}
	if(n->child == dead){
		/*this means this is the stored child and we know more exists since children > 1 */ 
		n->child = dead->rsib;
		dead->rsib->lsib = NULL;
	}
	else{
		dead->lsib->rsib = dead->rsib;
		if(dead->rsib != NULL){
			dead->rsib->lsib = dead->lsib;
		}
	}
	n->children--;
	zReleaseTBTreeNode(t,dead);	
	/* see if n is now redundant */
	if((n->children == 1) && (n->child->state == n->state)){
		zReleaseRedundantTBTreeNode(t,n);
	}
	/* we may have created a new cpoint here, if so find it */
	while(t->cpoint->children == 1){
		t->cpoint = t->cpoint->child;
		t->new_cpoint = 1;
		any_new_cpoint = t->cpoint->pos;
	}
}

/* this locks a node into the tree so that it is not removed if it becomes redundant */
void zTBTreeLockNode(zTBTreeNode* n){
	n->lock = 1;
}

/* release the lock on a node, removing it if it is redundant */
/*
static void zTBTreeUnLockNode(zTBTree* t, zTBTreeNode* n){
	n->lock = 0;
	if(n->parent != NULL && 
	   n->children == 1 && 
	   n->child->state == n->state){
		zReleaseRedundantTBTreeNode(t,n);
	}
}
*/

/* release all nodes which are decendants of *n in *t */ 
void zTBTreeClearSubTree(zTBTree* t, zTBTreeNode* n){
	zTBTreeNode* c = n->child;
	if(n->children == 0){
		while(c != NULL){
			zResetTBTreeHelp(t,c);
			c = c->rsib;
		}
		n->children = 0;
		n->child = NULL;
	}
}

void zTBTreeSetChild(zTBTreeNode* p,zTBTreeNode* c){
	/*remove from the old parent children */
	if(c->parent != NULL){
		if(c->parent->child == c){
			c->parent->child = c->rsib;
			/*if c is the stored child then it is the left most child (no lsib) */
		}
		else{
			/* otherwise it must have an lsib */
			c->lsib->rsib = c->rsib;
		}
		/* if right sibling exists, clear it */
		if(c->rsib != NULL){
			c->rsib->lsib = c->lsib;
		}
		c->parent->children--;
		if(c->parent->children < 0){
			zDie("zTBTreeLiveNode children count is off\n");
		}
	}
	c->lsib = NULL;
	c->rsib = NULL;
	/* make c p's leftmost child */
	if(p->child != NULL){
		c->rsib = p->child;
		p->child->lsib = c; 
	}
	p->child = c;
	c->parent = p;
	p->children++;
}

/* this copies everything under n1 in t1 to n2 in t2 
 n2 should have no children before this is called */
void zCopyTBTreeNode(zTBTree* t1, zTBTreeNode* n1, zTBTree* t2, zTBTreeNode* n2, int depth){
	zTBTreeNode* c1;
	zTBTreeNode* c2;
	zTBTreeNode* c3;
	if(t1->dead_nodes != NULL && t1->dead_nodes->alive == 1){
		zDie("");
	}
	if(n2->children != 0 || n1->alive == 0 || n2->alive == 0){
		zDie("");
	}
	n2->pos = n1->pos;
	n2->state = n1->state;
	n2->score = n1->score;
	n2->strand = n1->strand;
	n2->frame_data = n1->frame_data;
	n2->phase = n1->phase;
	n2->children = n1->children;
	c1 = n1->child;
	c3 = NULL;
	while(c1 != NULL){
		c2 = zGetTBTreeNode(t2);
		zCopyTBTreeNode(t1,c1,t2,c2,depth+1);
		c2->parent = n2;
		if(c3 == NULL){
			n2->child = c2;
			c2->lsib = NULL;
			c2->rsib = NULL;
		}
		else{
			c2->rsib = NULL;
			c2->lsib = c3;
			c3->rsib = c2;
		}
		c3 = c2;
		c1 = c1->rsib;
	}
}

/* this makes a copy of a zTBTree starting at the last cpoint */
void zCopyTBTreeFromCPoint(zTBTree* t1,zTBTree* t2,coor_t min_pos){	
	zTBTreeNode* n1 = t1->cpoint;
	zTBTreeNode* n2 = t2->root;
	while(n1 != t1->root && n1->pos > min_pos){
		n1 = n1->parent;
	}
	if(n1 != t1->root){
		n2 = zGetTBTreeNode(t2);
		zTBTreeSetChild(t2->root,n2);
		t2->root->pos = n1->parent->pos;
		t2->cpoint = n2;
	}
	zCopyTBTreeNode(t1,n1,t2,n2,0);
	/* update t2s cpiont */
	while(t2->cpoint->children == 1){
		t2->cpoint = t2->cpoint->child;
	}
}

/* copy everything from n1 down in t1 to node n2 in t2 */
void zAppendTBTreeFromNode(zTBTree* t1, zTBTreeNode* n1, zTBTree* t2, zTBTreeNode* n2){
	zCopyTBTreeNode(t1, n1, t2, n2,0);
	if(n2->parent != NULL && n2->children == 1 && n2->child->state == n2->state){
		zReleaseRedundantTBTreeNode(t2,n2);
	}
	else{
		zTBTreeNode* cp = n2;
		while(cp->children == 1){
			cp = cp->child;
		}
		t2->cpoint = cp;
	}
}
