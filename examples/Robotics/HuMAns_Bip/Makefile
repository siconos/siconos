
OBJS =	Contact.o \
	ContactJacobian.o \
	Inertia.o \
	JacobianQNLEffects.o \
	JacobianVNLEffects.o \
	NLEffects.o

all : dummy.so

%.o : %.c
	gcc -c -o $(@) -Wall -O3 $<

dummy.so : $(OBJS)
	gcc -shared $^ -o $@
	$(RM) dummy.so

