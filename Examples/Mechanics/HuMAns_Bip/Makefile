
OBJS =	Contact.o \
	ContactJacobian.o \
	ContactHessian.o \
	Inertia.o \
	JacobianQNLEffects1.o \
	JacobianQNLEffects2.o \
	JacobianQNLEffects3.o \
	JacobianQNLEffects4.o \
	JacobianQNLEffects5.o \
	JacobianQNLEffects6.o \
	JacobianQNLEffects7.o \
	JacobianQNLEffects.o \
	JacobianVNLEffects.o \
	NLEffects.o \
	Tags.o

all : libSiconosRobot.so

%.o : %.c
	gcc -c -o $(@) -Wall -O3 $<

libSiconosRobot.so : $(OBJS)
	gcc -shared $^ -o $@

