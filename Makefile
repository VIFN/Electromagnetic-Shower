rootlibs =$(shell root-config --libs)
rootflags =$(shell root-config --cflags)
HEADERINC = -I \headers
cc = $(shell root-config --cxx) -std=c++11 $(HEADERINC) $(rootflags) 

# corre propagator para 1 particula e cria cascata em 3d e histograma 
#(a aparecer 4 pads por predefinicao mas podem ser comentados se tal for conveniente)
# com "make e3 projecto" ou "make e5 projecto" é possivel alterar a energia da particula
# (default =1 GeV)
#com "make grafite projecto" corre o programa para a grafite. Outras opcoes sao iron e chumbo
#pode se combinar opcoes de energia e elemento

projecto: clean  particle.exe histogram.exe 3ddrawing.exe 
	@echo "\n" $< "\n"
	@./particle.exe
	@echo " made data \n"
	@./3ddrawing.exe
	@./histogram.exe

 ##### para fazer drawing 3d
full3d: clean 3ddrawing.exe particle.exe
	@echo "\n" $< "\n"
	@./particle.exe
	@echo " made data \n"
	@./3ddrawing.exe

3d: 3ddrawing.exe
	@echo "\n" $< "\n"
	@./$^


3ddrawing.exe: 3ddrawing.o  \libs/libFC.a
	@echo Linking $@ from $(notdir src/$^)
	$(cc) -o $@ $^  $(rootlibs)


### histogramas
fullhist: clean histogram.exe particle.exe
	@echo "\n" $< "\n"
	@./particle.exe
	@echo " made data \n"
	@./histogram.exe

histogram: histogram.exe
	@echo "\n" $< "\n"
	@./$^


histogram.exe: histogram.o  \libs/libFC.a 
	@echo Linking $@ from $(notdir src/$^)
	$(cc) -o $@ $^  $(rootlibs)


##criar dados (importante e funcional)
data: particle.exe
	@echo "\n" $< "\n"
	@./$^

particle.exe: mainParticle.o \libs/libPropagator.a \libs/libFormula.a \libs/libParticle.a \libs/libFCmatrix.a \libs/libVect.a \libs/libFC.a
	@echo Linking $@ from $(notdir src/$^)
	$(cc) -o $@ $^  $(rootlibs)

##histograma de nmax (apenas correr se tiver tempo para gastar ja que sao propagadas 2000 particulas por predefinicao)
histhmax.exe: nmaxhistogram.o \libs/libFC.a
	@echo Linking $@ from $(notdir src/$^)
	$(cc) -o $@ $^  $(rootlibs)

nmax.exe: nmaxtests.o \libs/libPropagator.a \libs/libFormula.a \libs/libParticle.a \libs/libFCmatrix.a \libs/libVect.a \libs/libFC.a
	@echo Linking $@ from $(notdir src/$^)
	$(cc) -o $@ $^  $(rootlibs)




particletest: particletest.exe
	@echo "\n" $< "\n"
	@./$^


particletest.exe: maintest.o \libs/libParticle.a  \libs/libVect.a \libs/libFC.a
	@echo Linking $@ from $(notdir src/$^)
	$(cc) -o $@ $^  $(rootlibs)

##### seria em D mas não funciona (razão desconhecida)
2d: 2ddrawing.exe
	@echo "\n" $< "\n"
	@./$^


2ddrawing.exe: 2ddrawing.o  \libs/libFC.a
	@echo Linking $@ from $(notdir src/$^)
	$(cc) -o $@ $^  $(rootlibs)

cFCgraphics.o: cFCgraphics.C
	@echo Compiling $(notdir src/$^)
	@$(cc) -c $< -o $@ 

%.o: %.cpp 
	@echo Compiling $(notdir src/$^)
	@$(cc) -c $< -o $@ 


debug: 
	$(eval cc += -DDEBUG)

wall:
	$(eval cc += -Wall)

iron: 
	$(eval cc += -DIRON)

chumbo: 
	$(eval cc += -DCHUMBO)

grafite: 
	$(eval cc += -DGRAFITE)

e3: 
	$(eval cc += -DE3)

e5: 
	$(eval cc += -DE5)


##### libs ###############################3

compilelib: \libs/libFormula.a \libs/libPropagator.a \libs/libParticle.a \libs/libFCmatrix.a \libs/libVect.a \libs/libFC.a 

\libs/libFCmatrix.a: FCmatrixFull.o FCmatrix.o
	@echo Compiling Library $(notdir src/$@)
	@ar rv $@ $^
	@ranlib $@


\libs/libFormula.a: Formula.o
	@echo Compiling Library $(notdir src/$@)
	@ar rv $@ $^
	@ranlib $@


\libs/libPropagator.a: Propagator.o
	@echo Compiling Library $(notdir src/$@)
	@ar rv $@ $^
	@ranlib $@

\libs/libVect.a: Vect.o
	@echo Compiling Library $(notdir src/$@)
	@ar rv $@ $^
	@ranlib $@

\libs/libFC.a: cFCgraphics.o
	@echo -n "Compiling" $@ "\n"
	@ar rv $@ $^

\libs/libParticle.a: Particle.o
	@echo Compiling Library $(notdir src/$@)
	@ar rv $@ $^
	@ranlib $@


# cleanups ##########################################

clean: cleanlocal cleanlib

cleanall: cleanlocal cleanlib cleanspecial cleanprints

cleanlocal:
	@echo Cleaning...
	@rm -f *.o
	@rm -f *.exe
	@rm -f *~


cleanspecial:
	@rm -f *.data

cleanprints:
	@rm -f *.png

cleanlib:
	@echo Cleaning...
	rm -f \libs/*.a




