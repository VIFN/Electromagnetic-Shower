Autores: Rafael Boto (84417), Vânia Nunes (85235)
Ao efetuar "make projecto", é corrido o propagator para 1 fotão de 1GeV e cria cascata em 3d e histograma (a aparecer 4 pads por predefinição mas podem ser comentados se tal for conveniente)
Com "make e3 projecto" ou "make e5 projecto" é possivel alterar a energia da particula (default =1 GeV)
Com "make grafite projecto" corre o programa para a grafite. Outras opções válidas sao iron e chumbo. Podem se combinar opcoes de energia e elemento
Codificou-se um programa nmaxtests.cpp e nmaxhistogram.cpp para obter histograma com o número máximo de partículas observadas em simultâneo num intervalo espacial, para determinada energia. 
Codificou-se um programa 2ddrawing.cpp que funcionaria de forma semelhante ao em 3D, mas utilizando os métodos do cFCgraphics.cpp propostos no FAQ deste projeto. No entanto, estes não funcionaram como previsto e considerámos não necessário a representação a 2 dimensões, já que se pode obter uma representação igual com o 3D em perspetiva.

Obrigado pela originalidade do projecto proposto já que conduziu-nos a aprender acerca da área da física envolvida, 
Rafael Boto & Vânia Nunes 
