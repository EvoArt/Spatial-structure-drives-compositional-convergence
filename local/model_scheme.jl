TikzPicture(L"""
                          \begin{scope}[
            yshift=-194,every node/.append style={
            yslant=0.5,xslant=-1},yslant=0.5,xslant=-1
            ]

        % opacity to prevent graphical interference
        \fill[white,fill opacity=1] (0,0) rectangle (6,6);        
         
        \fill[Askin] (4,4.5) rectangle (4.5,5.0); 
        \fill[Ahair] (2,4.5) rectangle (2.5,5.0); 
        %\draw[step=5mm, black] (0,0) grid (6,6);
        \draw[step=20mm, nobeadclr, thick] (0,0) grid (6,6); %defining
         \draw[black,very thick] (0,0) rectangle (6,6); 

		\fill[Ahair] (0,0) rectangle (0.5,0.5); 
        \fill[Ahair] (1.5,1.5) rectangle (2.0,2.0); 
        \fill[Adress2] (0,1.5) rectangle (0.5,2.0);       
		
        \fill[Ahair] (1.5,5) rectangle (2.0,5.5);         
		\fill[Rcollar] (2.5,3) rectangle (3,3.5);  
		
		\draw[thick] (0,4) rectangle (0.5,4.5); 
		\draw[thin, fill =Ahair] (-1.75,2.75) rectangle (-1.25,3.25); 		
		
		\draw[thick] (3,2) rectangle (2.5,1.5);
		\draw[thin, fill =Rcollar] (1.75,-0.75) rectangle (2.25,-0.25); 
		
		
		\draw[step=5mm, black] (0,0) grid (6,6); %defining grids
        %\draw[step=1mm, red!50,thin] (0.5,0.5) grid (2,2);  %Nested Grid
        \draw[black,very thick] (0,0) rectangle (6,6);%marking borders
		
		\draw[->,thick,red, bend right](0.25,4.25) to (-1.5,3);
		\draw[->,thick,red, bend left](2.75,1.75) to (2.0,-0.5);

    \end{scope}
   
             \begin{scope}[
            yshift=-97,every node/.append style={
            yslant=0.5,xslant=-1},yslant=0.5,xslant=-1
            ]

        % opacity to prevent graphical interference
        \fill[white,fill opacity=1] (0,0) rectangle (6,6);
        
		%\draw [thick,<->]  (0,6.3) -- node [yshift = 0.2cm, sloped] {\huge $d$} (2,6.3);        
        \fill[Ahair,fill opacity=0.3] (1,1) rectangle (1.5,1.5); 
        \fill[Ahair,fill opacity=0.3] (1.5,1) rectangle (2.0,1.5); 
        \fill[Rcollar,fill opacity=0.3] (3,2) rectangle (3.5,2.5); 
        \fill[Rcollar] (3,2) rectangle (2.5,1.5); 
        \fill[Adress2,fill opacity=0.3] (1.5,0.5) rectangle (2.0,1.0); 
        \fill[Askin] (4,4.5) rectangle (4.5,5.0); 
        \fill[Ahair,fill opacity=0.3] (0.5,5) rectangle (1.0,5.5); 
        \fill[Ahair] (2,4.5) rectangle (2.5,5.0); 
        \fill[Ahair,fill opacity=0.3] (2,4.5) rectangle (1.5,5.0); 
        %\draw[step=5mm, black] (0,0) grid (6,6);
       
		\fill[Ahair] (0,0) rectangle (0.5,0.5); 
        \fill[Ahair] (1.5,1.5) rectangle (2.0,2.0); 
        \fill[Adress2] (0,1.5) rectangle (0.5,2.0);       
        
		\fill[Ahair] (0,4) rectangle (0.5,4.5); 
        \fill[Ahair] (1.5,5) rectangle (2.0,5.5);         
				
		\fill[Rcollar] (2.5,3) rectangle (3,3.5);  		
		
		\draw[step=5mm, black] (0,0) grid (6,6); %defining grids
        %\draw[step=1mm, red!50,thin] (0.5,0.5) grid (2,2);  %Nested Grid
        \draw[black,very thick] (0,0) rectangle (6,6);%marking borders
        
         \draw[step=20mm, nobeadclr, thick] (0,0) grid (6,6); %defining
         \draw[black,very thick] (0,0) rectangle (6,6); 
         
         \draw[->,thick,nobeadclr, bend left](1.25,1.25) to (0.25,0.25);
		\draw[->,thick,nobeadclr](1.75,1.25) to (1.75,1.75);
		\draw[->,thick,nobeadclr, bend right](1.75,0.75) to (0.25,1.75);
		
		\draw[->,thick,nobeadclr, bend left](0.75,5.25) to (0.25,4.25);
		\draw[->,thick,nobeadclr](1.75,4.75) to (1.75,5.25);
		
		\draw[->,thick,nobeadclr, bend left](3.25,2.25) to (2.75,3.25);

    \end{scope}
           \begin{scope}[
            yshift=0,every node/.append style={
            yslant=0.5,xslant=-1},yslant=0.5,xslant=-1
            ]

        % opacity to prevent graphical interference
        \fill[white,fill opacity=1] (0,0) rectangle (6,6);
        %\draw[step=5mm, black] (0,0) grid (6,6); %defining grids
        %\draw[step=1mm, red!50,thin] (0.5,0.5) grid (2,2);  %Nested Grid
        \draw[black,very thick] (0,0) rectangle (6,6);%marking borders
        \fill[Ahair] (1,1) rectangle (1.5,1.5); 
        \fill[Ahair] (1.5,1) rectangle (2.0,1.5); 
        \fill[Rcollar] (3,2) rectangle (3.5,2.5); 
        \fill[Rcollar] (3,2) rectangle (2.5,1.5); 
        \fill[Adress2] (1.5,0.5) rectangle (2.0,1.0); 
        \fill[Askin] (4,4.5) rectangle (4.5,5.0); 
        \fill[Ahair] (0.5,5) rectangle (1.0,5.5); 
        \fill[Ahair] (2,4.5) rectangle (2.5,5.0); 
        \fill[Ahair] (2,4.5) rectangle (1.5,5.0); 
        \draw[step=5mm, black] (0,0) grid (6,6); %defining grids

    \end{scope}
       \begin{scope}[
            yshift=97,every node/.append style={
            yslant=0.5,xslant=-1},yslant=0.5,xslant=-1
            ]

        % opacity to prevent graphical interference
        \fill[white,fill opacity=1] (0,0) rectangle (6,6);
        \draw[step=5mm, black] (0,0) grid (6,6); %defining grids
        %\draw[step=1mm, red!50,thin] (0.5,0.5) grid (2,2);  %Nested Grid
        \draw[black,very thick] (0,0) rectangle (6,6);%marking borders
        %\draw[black,very thick,<-> ] (0,6.3) to (6,6.3) node [center, sloped] {\huge $S$} ;%marking borders
        %\draw [thick,<->]  (0,6.3) -- node [yshift = 0.2cm, sloped] {$S$} (6,6.3);
        \fill[Ahair] (1,1) rectangle (1.5,1.5); 
        \fill[Rcollar] (3,2) rectangle (3.5,2.5); 
        \fill[Adress2] (1.5,0.5) rectangle (2.0,1.0); 
        \fill[Askin] (4,4.5) rectangle (4.5,5.0); 
        \fill[Ahair] (0.5,5) rectangle (1.0,5.5); 
        \fill[Ahair] (2,4.5) rectangle (2.5,5.0); 
  
    \end{scope}
    
                \draw[-latex,thick](7,5.8)node[right,align=left]{\huge \textbf{Evaluate} fitness of each individual.}
        to[out=180,in=90] (5,5.8);
        
     		            \draw[-latex,thick](7,2.4)node[right,align=left]{\huge \textbf{Births} take place on neighbouring\\ \huge tiles.}
        to[out=180,in=90] (5,2.4);
  
            \draw[-latex,thick](7,-1.0)node[right,align=left]{\huge Individuals \textbf{disperse} to random tiles.}
        to[out=180,in=90] (5,-1);
        
                   \draw[-latex,thick](7,-4.4)node[right,align=left]{\huge Individuals are \textbf{removed} at random\\ \huge until $\text{density} = \text{limit}$.}
        to[out=180,in=90] (5,-4.4);
""", preamble = "\\usepackage{amsmath}
\\definecolor{nobeadclr}{RGB}{0,47,74}
\\definecolor{liquidclr}{RGB}{214,42,41}
\\definecolor{beadclr}{RGB}{253,192,76}
\\definecolor{Ahair}{RGB}{244,237,70}
\\definecolor{Adress1}{RGB}{180,242,255}
\\definecolor{Adress2}{RGB}{92,225,244}
\\definecolor{Askin}{RGB}{255,234,202}
\\definecolor{Rcollar}{RGB}{255,65,65}
\\usepackage[scaled]{helvet}
\\usepackage[T1]{fontenc}
\\renewcommand\\familydefault{\\sfdefault}
"
)