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

#ifndef __gmatrixH
#define __gmatrixH

#include <string> 
#include <vector>
#include <cstdlib>
#include "ebasics.h"

using namespace std;

struct mItem{
	int M;   // M-matrix item
	int Ix;  // Ix - gap matrix
	int Iy;  // Iy - gap matrix
	int GapH; // gap length ( Ix )
	int GapV; // gap length ( Iy )
};

typedef vector< vector< mItem> > mMatrix; // matrix object

struct mAlignment{
	string a1;
	string a2;
	string mL;
	int score;
	int a1_start;
	int a2_start;
	int a1_end;
	int a2_end;
	
	int positives;
	int identities;
	int gaps; 
	int a1_alignedP;
	int a2_alignedP;
};

class generalmatrix{
	public:
		generalmatrix() {}
		~generalmatrix() {}
		
		// matrix 
		void clearmatrix();
		void buildmatrix();
		
		// accession to matrix
		
		void setM( int i, int j, int M ) { Matrix[i][j].M = M; }
		void setIx( int i, int j, int Ix ) { Matrix[i][j].Ix = Ix; }
		void setIy( int i, int j, int Iy ) { Matrix[i][j].Iy = Iy; }
		void setHgap( int i, int j, int gap ) { Matrix[i][j].GapH = gap; }
		void setVgap( int i, int j, int gap ) { Matrix[i][j].GapV = gap; }
		
		int getM( int i, int j ) { return Matrix[i][j].M; }
		int getIx( int i, int j ) { return Matrix[i][j].Ix; }
		int getIy( int i, int j ) { return Matrix[i][j].Iy; }
		
		int iMax() { return iCols - 1; }
		int jMax() { return iRows - 1; }
		
		int getV( int i, int j ) { return max( max( getIx( i, j ), getIy( i, j ) ), getM( i, j ) ); }
		
		int getHgap( int i, int j ) { return Matrix[i][j].GapH; }
		int getVgap( int i, int j ) { return Matrix[i][j].GapV; }
		
		int getT( int i, int j );

		// alignment issues
		void initialfill();
		void matrixfill();
		void traceback();
		
		virtual int score( int i, int j ) { return 0; }
		virtual char mType( int i, int j ) { return ' '; }
		// matrix accession
		mAlignment getalignment() { return alignment; }
		
		void setPenalty( int gO, int gE ) { gapOpen = gO; gapExtention = gE; }
		void setSequences( generalSequence S1, generalSequence S2);
		
		void doalignment();
		void setmodus( bool lonly, bool coverhang, bool showskip);
		
		// scores for gap/identity/positives
		virtual void calcScores() {} 
		virtual void display_specific() {} 
		
		void dumpAlignment( int colW );
		void dumpAlignment_COMPUTERS(); // new 13-10-2006
		
		void calcselfsims();
		
		// outsiders acession 17-8-2006
		
		int nCols() { return iCols; }
		int nRows() { return iRows; }
		int vOpen() { return gapOpen; }
		int vExtention() { return gapExtention; }
		bool mode_localOnly() { return localOnly; }
		bool mode_chargeOverHang() { return chargeOverHang; } 
		bool mode_showSkip() { return showSkip; }
		
		int locali() { return local_i; }
		int localj() { return local_j; }
		int overhangi() { return overhang_i; }
		int overhangj() { return overhang_j; }
		
		generalSequence sSequence1() { return Sequence1; }
		generalSequence sSequence2() { return Sequence2; }
		
		string Title() { return myTitle; } // 13-10-2006
		void setTitle( string title ) { myTitle = title; } // 13-10-2006
		
	protected:
		mMatrix Matrix;
		mAlignment alignment;
		int iCols;
		int iRows;
		
		generalSequence Sequence1;
		generalSequence Sequence2;
		
		int gapOpen;
		int gapExtention; 
		
		bool localOnly;
		bool chargeOverHang;
		bool showSkip;
		
		int local_j;
		int local_i;
		
		int overhang_i;
		int overhang_j;
		
		float selfSc_1; 
		float selfSc_2;
		
		float sRatio_1;
		float sRatio_2; 
		
		string myTitle; // 13-10-2006
		
		

};
	
	
#endif 