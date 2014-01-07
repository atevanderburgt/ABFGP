/***************************************************************************
 *   Copyright (C) 2006 by edouard severing   *
 *   edouard@localhost   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef __alignprotsH
#define __alignprotsH

#include "gmatrix.h"

typedef int proteinMatrix[26][26];

class alignproteins: public generalmatrix{

	public:
		alignproteins(){}
		~alignproteins(){}
		
		int score( int i, int j );
		char mType( int i, int j ); 
		void loadMatrix( proteinMatrix pM );
		void loadInHouseMatrix( std::string mName); 
		
		void calcScores(); 
		void display_specific();

	protected:
	
		proteinMatrix pMatrix;
		std::string MatrixName;
		
};


extern proteinMatrix blosum100;
extern proteinMatrix blosum30;
extern proteinMatrix blosum35;
extern proteinMatrix blosum40;
extern proteinMatrix blosum45;
extern proteinMatrix blosum50;
extern proteinMatrix blosum55;
extern proteinMatrix blosum60;
extern proteinMatrix blosum62;
extern proteinMatrix blosum65;
extern proteinMatrix blosum70;
extern proteinMatrix blosum75;
extern proteinMatrix blosum80;
extern proteinMatrix blosum85;
extern proteinMatrix blosum90;
extern proteinMatrix blosum95;
extern proteinMatrix pam120;
extern proteinMatrix pam180;
extern proteinMatrix pam250;
extern proteinMatrix pam30;
extern proteinMatrix pam300;
extern proteinMatrix pam60;
extern proteinMatrix pam90;
extern proteinMatrix blosum50;

#endif
