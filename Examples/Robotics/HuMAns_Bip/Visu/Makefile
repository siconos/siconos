include MakefileSys

LIB = libLagrangianModel$(LIBEXT)



BIPOBJS = Tags.o
       

OBJS=$(BIPOBJS) ScilabInterface.o

all : $(OBJPATH)/$(LIB) 

$(OBJPATH)/%.o: %.c
	$(CC) -c -o $(@) $(CFLAGS) $<

$(OBJPATH)/$(LIB) : $(addprefix $(OBJPATH)/, $(OBJS))
	$(LD) $(LDFLAGS) $^ -o $@

clean :
	$(RM) $(OBJPATH)/*.o $(OBJPATH)/$(LIB) 

SCIFILES=SomeDefinitions.sci


