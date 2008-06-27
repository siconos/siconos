include MakefileSys

LIB = libLagrangianModel$(LIBEXT)



BIPOBJS = Contact.o \
	ContactHessian.o \
	ContactJacobian.o \
	Inertia.o \
	NLEffects.o \
	JacobianNLEffects.o \
	JacobianVelocityNLEffects.o \
	Tags.o 

OBJS=$(BIPOBJS) ScilabInterface.o

all : $(OBJPATH)/$(LIB) 

$(OBJPATH)/%.o: %.c
	$(CC) -c -o $(@) $(CFLAGS) $<

$(OBJPATH)/$(LIB) : $(addprefix $(OBJPATH)/, $(OBJS))
	$(LD) $(LDFLAGS) $^ -o $@


clean :
	$(RM) $(OBJPATH)/*.o $(OBJPATH)/$(LIB) 

SCITOC=scitoc
SCIFILES=SomeDefinitions.sci

CFILES=$(patsubst %.sci,%.scilab,$(SCIFILES))

doc: $(CFILES)


$(CFILES): %.scilab : %.sci
	$(SCITOC) $< >$(@)

clean_doc:	
	$(RM) *.scilab

bip: all
	$(CP) $(OBJPATH)/$(BIPLIB) ../../Tools/BipExperiments/Export/$(OBJPATH)
	$(CP) LagrangianModel.h ../../Tools/BipExperiments/Export/include
