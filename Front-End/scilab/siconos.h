
extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicLoadModel(int *ret, char ModelXmlFile[]);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicInitStrategy();

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicTimeGetH(double *H);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicTimeGetN(int *N);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicTimeGetK(int *K);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void  sicSTNextStep(int *ret);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicSTComputeFreeState(int *ret);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicSTcomputePb(int *ret);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicSTupdateState(int *ret);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicDebug(int *ret);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicModelgetQ(double *value, int *index);

