#ifndef SUPPORT_CCRANDSTATE_CC_H
#define SUPPORT_CCRANDSTATE_CC_H

int blub();
#define CC_PRANDMAX 1000000007

typedef struct CCrandstate {
			int a;
			int b;
			int arr[55];
		} CCrandstate;

int
	CCutil_lprand (CCrandstate *r);

void
   CCutil_sprand (int seed, CCrandstate *r);

#endif
