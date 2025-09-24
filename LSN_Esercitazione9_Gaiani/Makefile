CC = g++
CFLAGS = -std=c++14 -Wall -Wextra
LDLIBS = -larmadillo

TARGET = main.exe
OBJS = main.o random.o

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDLIBS)

main.o: main.cpp population.h individual.h city.h random.h
	$(CC) $(CFLAGS) -c main.cpp

random.o: random.cpp random.h
	$(CC) $(CFLAGS) -c random.cpp

clean:
	rm -f *.o $(TARGET) seed.out
