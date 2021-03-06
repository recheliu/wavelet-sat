#pragma once

#if		!defined(WITH_DYNAMIC_ARRAY_ALLOCATION)
	#define WITH_DYNAMIC_ARRAY_ALLOCATION	1
#endif	// #if	!defined(WITH_DYNAMIC_ARRAY_ALLOCATION)

#if		WITH_DYNAMIC_ARRAY_ALLOCATION	
	#if		!defined(WITH_STREAMING)
		#define WITH_STREAMING	1
	#endif	// #if	!defined(WITH_STREAMING)
#endif	// #if		WITH_DYNAMIC_ARRAY_ALLOCATION	

#include <map>	
#if defined (WIN32)
	#include <psapi.h>	
	#pragma comment (lib, "psapi.lib")
#else	// #if defined (WIN32)
	#include <sys/time.h>
	#include <sys/resource.h>
#endif	// #if defined (WIN32)

#include <vector>
using namespace std;
#include <math.h>

#include "HeaderBase.h"	
#include "DWT.h"
#include "liblog.h"	

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The base class of sep. DWT 
	/*!
	In order to apply wavelet transform for SATs or Integral Histograms, please follow the procedure below:
	1. Setup the data size by _SetDimLengths().
	2. Setup the #SATs by _AllocateBins().
	3. Call _Update() for each data point to incrementally compute wavelet.
	4. Call _Finalize() to finish the computation.

	To query the entry for a certain location:
	1. Call _GetAllSums() to get all SAT values for the given location.
	*/
	class CSepDWTHeader:
		virtual public CHeaderBase
	{
protected:	
		//! Total #coefficients
		size_t uNrOfCoefs;

		//! #Coefs per dim.
		/*!
		c[0], ..., c[d], ... c[D - 1]
		*/
		vector<size_t> vuCoefLengths;

		//! The weight sum for each bin
		vector<double> vdBinWeights;

		//! The maximun number of wavelet levels per dim, which are specified by the user.
		/*!
		m[0], ... m[d], ... m[D - 1]
		*/
		vector<size_t> vuDimMaxLevels;
		
		//! The total number of wavelet basis from all dimensions.
		/*!
		m[0] + m[d] + ... m[D - 1]
		*/
		size_t uNrOfWaveletsFromAllDims;

		//! #levels per dim
		/*!
		l[0] = log2(n[0]) + 1, ... l[d] = log2(n[d]) + 1, ... l[D - 1] = log2(n[D - 1]) + 1
		*/
		vector<size_t> vuDimLevels;

		//! Product of influencing numbers of coefficients from all dimensions
		/*! 
		#coefficients to update ...
		C = prod l[d]
		*/
		size_t uNrOfUpdatingCoefs;

		vector< size_t >			vuMapLocalToGlobal;	//!< Map 1D local index to its global index
		vector< vector<size_t> >	vvuGlobalBase;		//!< D-dim global base foe each Wavelet
		vector< size_t >			vuGlobalBase;		//!< 1-dim global base foe each Wavelet
		vector< vector<size_t> >	vvuLocalLengths;	//!< vector of D-dim tuples as the coefficient lengths
		vector< size_t >			vuNrsOfLocalCoefs;	//!< vector of #coefficients per wavele

		//! vector of D-dim tuples as the lengths fo the wavelet functions.
		vector< vector<size_t> >	vvuWaveletLengths;	

		//! vector of D-dim tuples as the coefficient lengths within the data.
		vector< vector<size_t> >	vvuLocalDataLengths;

		//! The base of all scanlines for all dimensions
		vector< vector<size_t> > vvuSliceScanlineBase;

		vector< vector<bool> > vvbSliceScanlineValid;	

		//! A lookup table to map the coefficient to its levels per dim.
		/*! size: D x C: 
		Map each coefficients c, c = 0 ... C - 1 = prod l[d] - 1, to the the corresponding D levels
		*/
		vector<size_t> vuCoefDim2Level;

		//! A lookup table to map the coefficient to its 1D index in the L x D wavelets
		/*! size: D x C: 
		Map each coefficients c, c = 0 ... C - 1 = prod l[d] - 1, to the the corresponding index 
		*/
		vector<size_t> vuCoefDim2Wavelet;

		//! The look up table per dimension to quickly conver dim. subscript to the coefficient indices
		/*!
		For each dim d, 
			the size is l[d] x n[d]: map the spatial location i, i = 0 ... n[d], to its l[d] coefficient locations
		*/
		vector< vector<size_t> > vvuSubLevel2Coef;

		//! Max. count per coefficients
		vector<size_t> vuMaxCounts;	

		//! The denomator for Wavelet basis
		/*!
		W = sqrt(prod n[0], ... n[d], ... n[D])
		*/
		double dWaveletDenomiator;

		//! The mininal value of the projected wavelet coefficients
		double dWaveletThreshold;	
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0100,
			PARAMETER_END
		};

		struct CDimLevel {
			size_t uWaveletLength;
			size_t uCoefLength;
			size_t uDataLength;
		};

		vector< vector<CDimLevel> > vvcDimLevels;

		const vector<size_t>& 
		VGetDimLevels
		(
			void *_Reserved = NULL
		) const 
		{	
			return vuDimLevels;
		}

		const vector<size_t>& 
		VGetCoefLengths
		(
			void *_Reserved = NULL
		) const 
		{	
			return this->vuCoefLengths;
		}

		//! Retunr the smallest value. 
		/*!
		This can be used to filter too small value caused by numerical error.
		*/
		virtual	
		double 
		DGetThreshold
		(
			void *_Reserved = NULL
		)
		{
			return dWaveletThreshold;
		}

		/*!
		Convert the (global) index of a wavelet basis to its corresponding subscripts and levels. 
		Also, for the wavelet function for that level, return the global base, local coef. lengths, 
		and the offset in this level.
		*/
		virtual
		void
		_ConvertIndexToLevels
		(
			size_t uIndex,
			vector<size_t>& vuLevel,
			vector<size_t>& vuLocalCoefSub,
			vector<size_t>& vuGlobalCoefBase,
			vector<size_t>& vuLocalCoefLengths,
			void *_Reserved = NULL
		)
		{
			vector<size_t> vuSub;
			_ConvertIndexToSub(uIndex, vuSub, vuCoefLengths);

			if( UGetNrOfDims() != vuLevel.size() )				vuLevel.resize(UGetNrOfDims());
			if( UGetNrOfDims() != vuGlobalCoefBase.size() )		vuGlobalCoefBase.resize(UGetNrOfDims());
			if( UGetNrOfDims() != vuLocalCoefLengths.size() )	vuLocalCoefLengths.resize(UGetNrOfDims());
			if( UGetNrOfDims() != vuLocalCoefSub.size() )		vuLocalCoefSub.resize(UGetNrOfDims());
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				vuLevel[d] = (!vuSub[d])?0:(size_t)(1 + floor(log( (double)vuSub[d]) / M_LN2));
			}
			_ConvertWaveletSubToLevels(vuLevel, vuGlobalCoefBase, vuLocalCoefLengths);
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				vuLocalCoefSub[d] = vuSub[d] - vuGlobalCoefBase[d];
			}
		}

		virtual
		void
		_ConvertWaveletToLevels
		(
			size_t	uWavelet,
			vector<size_t>& vuGlobalCoefBase,
			vector<size_t>& vuLocalCoefLengths,
			size_t& uNrOfLocalCoefs,
			void *_Reserved = NULL
		)
		{
			if( UGetNrOfDims() != vuGlobalCoefBase.size() )
				vuGlobalCoefBase.resize(UGetNrOfDims());
			if( UGetNrOfDims() != vuLocalCoefLengths.size() )
				vuLocalCoefLengths.resize(UGetNrOfDims());
			copy(this->vvuGlobalBase[uWavelet].begin(), this->vvuGlobalBase[uWavelet].end(), vuGlobalCoefBase.begin());
			copy(this->vvuLocalLengths[uWavelet].begin(), this->vvuLocalLengths[uWavelet].end(), vuLocalCoefLengths.begin());
			uNrOfLocalCoefs = vuNrsOfLocalCoefs[uWavelet];
		}

		virtual
		void
		_ConvertWaveletSubToLevels
		(
			const vector<size_t>&	vuWaveletSub,
			vector<size_t>& vuGlobalCoefBase,
			vector<size_t>& vuLocalCoefLengths,
			void *_Reserved = NULL
		)
		{
			size_t uNrOfLocalCoefs;
			_ConvertWaveletToLevels(UConvertSubToIndex(vuWaveletSub, vuDimLevels), vuGlobalCoefBase, vuLocalCoefLengths, uNrOfLocalCoefs);
		}

		//////////////////////////////////////////////////////////////////////////
		const
		size_t 
		UGetDimLevelWaveletLength
		(
			const size_t uDim,
			const size_t uLevel,
			void *_Reserved = NULL
		) 
		const	
		{
			ASSERT_OR_LOG(uDim < UGetNrOfDims(), "");
			ASSERT_OR_LOG(uLevel < vvcDimLevels[uDim].size(), "");
			return vvcDimLevels[uDim][uLevel].uWaveletLength;
		}

		size_t UGetDimLevelCoefLength
		(
			const size_t uDim,
			const size_t uLevel,
			void *_Reserved = NULL
		)
		{
			ASSERT_OR_LOG(uDim < UGetNrOfDims(), "");
			ASSERT_OR_LOG(uLevel < vvcDimLevels[uDim].size(), "");
			return vvcDimLevels[uDim][uLevel].uCoefLength;
		}

		size_t UGetDimLevelDataLength
		(
			const size_t uDim,
			const size_t uLevel,
			void *_Reserved = NULL
		)
		{
			ASSERT_OR_LOG(uDim < UGetNrOfDims(), "");
			ASSERT_OR_LOG(uLevel < vvcDimLevels[uDim].size(), "");
			return vvcDimLevels[uDim][uLevel].uDataLength;
		}

		//////////////////////////////////////////////////////////////////////////
		const vector<size_t>&
		VGetGlobalCoefBase
		(
			const size_t	uWavelet,
			void *_Reserved = NULL
		) const 
		{
			return this->vvuGlobalBase[uWavelet];
		}

		const vector<size_t>&
		VGetGlobalCoefBase
		(
			const vector<size_t>&	vuWaveletSub,
			void *_Reserved = NULL
		) const 
		{
			return this->vvuGlobalBase[UConvertSubToIndex(vuWaveletSub, vuDimLevels)];
		}

		//////////////////////////////////////////////////////////////////////////
		const vector<size_t>&
		VGetLocalLengths
		(
			const size_t	uWavelet,
			void *_Reserved = NULL
		) const 
		{
			return this->vvuLocalLengths[uWavelet];
		}

		const vector<size_t>&
		VGetLocalLengths
		(
			const vector<size_t>&	vuWaveletSub,
			void *_Reserved = NULL
		) const 
		{
			return this->vvuLocalLengths[UConvertSubToIndex(vuWaveletSub, vuDimLevels)];
		}

		//////////////////////////////////////////////////////////////////////////
		const vector<size_t>&
		VGetWaveletLengths
		(
			const size_t	uWavelet,
			void *_Reserved = NULL
		) const 
		{
			return this->vvuWaveletLengths[uWavelet];
		}

		const vector<size_t>&
		VGetWaveletLengths
		(
			const vector<size_t>&	vuWaveletSub,
			void *_Reserved = NULL
		) const 
		{
			return this->vvuWaveletLengths[UConvertSubToIndex(vuWaveletSub, vuDimLevels)];
		}

		//////////////////////////////////////////////////////////////////////////
		const vector<size_t>&
		VGetLocalDataLengths
		(
			const size_t	uWavelet,
			void *_Reserved = NULL
		) const 
		{
			return this->vvuLocalDataLengths[uWavelet];
		}

		const vector<size_t>&
		VGetLocalDataLengths
		(
			const vector<size_t>&	vuWaveletSub,
			void *_Reserved = NULL
		) const 
		{
			return this->vvuLocalDataLengths[UConvertSubToIndex(vuWaveletSub, vuDimLevels)];
		}

		virtual
		void 
		_GetForwardWavelet(
			const vector<size_t>& vuPos, 
			vector<size_t>& vuSubs, 
			vector<long>& vlWavelets, 
			bool bIsNorminatorOnly,
			void* _Reserved = NULL)
		{
			vuSubs.resize(uNrOfWaveletsFromAllDims);
			vlWavelets.resize(uNrOfWaveletsFromAllDims);
			for(size_t p = 0, d = 0; d < UGetNrOfDims(); d++)
			{
				size_t uDimMaxLevel = vuDimMaxLevels[d];
				size_t uPos = vuPos[d];
				size_t uMaxWin = 1 << vuDimLevels[d];
				for(size_t 	
					l = 0, w = uMaxWin;
					l < uDimMaxLevel; 
					l++, p++, w >>= 1)
				{
					long lWavelet;
					// Given the currenet level l and subscript uPos, compute the sum of the portion in wavelet after uPos
					size_t uPosInWavelet = uPos % w;
					if( 0 == l )
						lWavelet = (long)w / 2 - (long)uPosInWavelet;
					else
					{
						if( uPosInWavelet < w / 2)
							lWavelet = (long)uPosInWavelet;
						else
							lWavelet = (long)(w - uPosInWavelet);
						lWavelet *= -1;		
					}
					vlWavelets[p] = lWavelet;
					vuSubs[p] = uPos / w;
				}
			}
		}

		virtual
		void
		_GetBackwardWavelet(
			const vector<size_t>& vuPos, 
			vector<size_t>& vuSubs, 
			vector<double>& vdWaveletBasis, 
			bool bIsNorminatorOnly,
			void* _Reserved = NULL)
		{
			vdWaveletBasis.resize( this->uNrOfWaveletsFromAllDims );
			vuSubs.resize( this->uNrOfWaveletsFromAllDims );
			// for each dimenion d, based on the posistion, store the corresponding l[d] wavelet basis value
			for(size_t p = 0, d = 0; d < UGetNrOfDims(); d++)
			{
				size_t uPos = vuPos[d];
				size_t uDimMaxLevel = vuDimMaxLevels[d];
				size_t uMaxWin = 1 << vuDimLevels[d];
				for(size_t 	
					l = 0, w = uMaxWin; 
					l < uDimMaxLevel; 
					l++, p++, w >>= 1)
				{
					// Decide the wavelet size based on the current level, and then the wavelet basis values based on the position within this wavelet
					double dWaveletBasis = 0.0;
					
					// assume Haar wavelet for now
					// Given the currenet level l and subscript uPos, compute the corresponding wavelet value
					size_t uPosInWavelet = uPos % w;				
					if( uPosInWavelet < w / 2)
						dWaveletBasis = 1.0;
					else
						dWaveletBasis = -1.0;

					if( l >= 2 )
						dWaveletBasis *= sqrt( (double)(1 << (l - 1)) );
					vdWaveletBasis[p] = dWaveletBasis;

					vuSubs[p] = uPos / w;
				}
			}
		}

		virtual	
		void 
		_SetMaxLevels
		(
			vector<size_t>& vuDimMaxLevels,
			void *_Reserved = NULL
		)
		{
			this->vuDimMaxLevels.clear();
			uNrOfCoefs = 1;
			this->uNrOfUpdatingCoefs = 1;
			this->uNrOfWaveletsFromAllDims = 0;	
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				size_t uMaxLevel = vuDimMaxLevels[d];
				size_t uCoefLength;	// total number of coefficient for this dimension
				if( !uMaxLevel )
					uMaxLevel = vuDimLevels[d];
				else
					uMaxLevel = min(uMaxLevel, vuDimLevels[d]);
				this->uNrOfUpdatingCoefs *= uMaxLevel;
				uCoefLength = (size_t)1 << (uMaxLevel - 1);
				vuCoefLengths[d] = uCoefLength;	
				this->vuDimMaxLevels.push_back(uMaxLevel);
				uNrOfCoefs *= uCoefLength;		// update the total number of coefficient to store for all dimension

				uNrOfWaveletsFromAllDims += uMaxLevel;	
			}

			this->vuCoefDim2Level.resize(uNrOfUpdatingCoefs * UGetNrOfDims());
			this->vuCoefDim2Wavelet.resize(uNrOfUpdatingCoefs * UGetNrOfDims());
			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
				for(size_t	uCoef = c, uBase = 0,	d = 0; 
						d < UGetNrOfDims(); 
						uBase += this->vuDimMaxLevels[d], d++)
				{
					size_t uMaxLevel = this->vuDimMaxLevels[d];
					size_t uLevel = uCoef % uMaxLevel; 
					uCoef /= uMaxLevel; 
					this->vuCoefDim2Level[c * UGetNrOfDims() + d] = uLevel;
					this->vuCoefDim2Wavelet[c * UGetNrOfDims() + d] = uBase + uLevel;
 				}

			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				size_t uMaxCountPerCoef = 1;	// max # coefficnet
				for(size_t d = 0; d < UGetNrOfDims(); d++, p++)
				{
					size_t uLevel = vuCoefDim2Level[p];
					uMaxCountPerCoef *= (size_t)1 << (( 0 == uLevel )?(vuDimLevels[d] - 1 - uLevel):(vuDimLevels[d] - uLevel));
				}
				vuMaxCounts.push_back(uMaxCountPerCoef);
			}

			vvcDimLevels.resize(UGetNrOfDims());
			for(size_t d = 0; d < UGetNrOfDims(); d++) 
			{
				vvcDimLevels[d].resize(vuDimLevels[d]);
				for(size_t l = 0; l < vuDimLevels[d]; l++) 
				{
					size_t uWaveletLength = (!l)?( 1<<(vuDimLevels[d] - 1) ):(1<<(vuDimLevels[d] - l)); 
					vvcDimLevels[d][l].uWaveletLength = uWaveletLength;
					vvcDimLevels[d][l].uDataLength = (size_t)ceil((double)vuDimLengths[d]/(double)uWaveletLength);
					vvcDimLevels[d][l].uCoefLength = (!l)?1:(1<<(l-1));
				}
			}

			vvuGlobalBase.resize(uNrOfUpdatingCoefs);
			vuGlobalBase.resize(uNrOfUpdatingCoefs);
			vvuLocalLengths.resize(uNrOfUpdatingCoefs);
			vuNrsOfLocalCoefs.resize(uNrOfUpdatingCoefs);
			vuMapLocalToGlobal.resize(uNrOfCoefs);
			vvuWaveletLengths.resize(uNrOfUpdatingCoefs);
			vvuLocalDataLengths.resize(uNrOfUpdatingCoefs);

			vector<size_t> vuLevel;
			for(size_t i = 0, w = 0; w < uNrOfUpdatingCoefs; w++)
			{
				_ConvertIndexToSub(w, vuLevel, vuDimLevels);
				vvuGlobalBase[w].resize(vuLevel.size());
				vvuLocalLengths[w].resize(vuLevel.size());
				vvuWaveletLengths[w].resize(vuLevel.size());
				vvuLocalDataLengths[w].resize(vuLevel.size());
				size_t uLength = 1;
				for(size_t d = 0; d < vuLevel.size(); d++)
				{
					size_t uLevel = vuLevel[d];
					vvuGlobalBase[w][d] = (!uLevel)?0:(1 << (uLevel - 1));
					vvuLocalLengths[w][d] = (!uLevel)?1:(1 << (uLevel - 1));
					vvuWaveletLengths[w][d] = UGetDimLevelWaveletLength(d, uLevel);
					vvuLocalDataLengths[w][d] = UGetDimLevelDataLength(d, uLevel);
					uLength *= vvuLocalLengths[w][d];
				}
				vuNrsOfLocalCoefs[w] = uLength;

				vector<size_t> vuLocalCoef;		vuLocalCoef.resize(UGetNrOfDims());
				vector<size_t> vuGlobalCoef;	vuGlobalCoef.resize(UGetNrOfDims());
				for(size_t lc = 0; lc < uLength; lc++, i++)
				{
					_ConvertIndexToSub(lc, vuLocalCoef, vvuLocalLengths[w]);
					for(size_t d = 0; d < UGetNrOfDims(); d++)
						vuGlobalCoef[d] = vvuGlobalBase[w][d] + vuLocalCoef[d];
					vuMapLocalToGlobal[i] = UConvertSubToIndex(vuGlobalCoef, vuCoefLengths);
				}
			}

			vector<size_t> vuLocalCoefSub, vuGlobalCoefBase, vuLocalCoefLengths;
			vvbSliceScanlineValid.resize(UGetNrOfDims());
			size_t uNrOfInvalids = 0;

			vvuSliceScanlineBase.resize(UGetNrOfDims());
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				vector<size_t> vuSliceCoefLengths = vuCoefLengths;
				vuSliceCoefLengths[d] = 1;
				size_t uNrOfScanlines = uNrOfCoefs / vuCoefLengths[d];
				vvuSliceScanlineBase[d].resize(uNrOfScanlines);
				vvbSliceScanlineValid[d].resize(uNrOfScanlines);
				for(size_t s = 0; s < uNrOfScanlines; s++)
				{
					vector<size_t> vuScanlineBase;
					_ConvertIndexToSub(s, vuScanlineBase, vuSliceCoefLengths);
					vvuSliceScanlineBase[d][s] = UConvertSubToIndex(vuScanlineBase, vuCoefLengths);
					this->_ConvertIndexToLevels(
						vvuSliceScanlineBase[d][s],
						vuLevel,
						vuLocalCoefSub,
						vuGlobalCoefBase,
						vuLocalCoefLengths);
					bool bIsValid = true;
					for(size_t d2 = 0; d2 < UGetNrOfDims(); d2++)
						if( vuLevel[d2] > 0 )
						{
							size_t uWaveletLength = vuCoefLengths[d2] / ((size_t)1 << (vuLevel[d2] - 1));
							// uWaveletLength /= 2;	// TMP
							if( vuLocalCoefSub[d2] >= (size_t)ceil( (double)vuDimLengths[d2]/(double)uWaveletLength ) )
							{
								bIsValid = false;
								uNrOfInvalids++;
								break;
							}
						}

					vvbSliceScanlineValid[d][s] = bIsValid;	// (!d)?bIsValid:true;
				}
			}
		}

		//! Set up the data dimensions
		virtual
		void
		_Set
		(
			const vector<size_t>& vuDimLengths,
			size_t uNrOfBins,
			void *_Reserved = NULL
		)
		{
			CHeaderBase::_Set(vuDimLengths, uNrOfBins);

			vdBinWeights.resize(uNrOfBins);

			this->vuCoefLengths.clear();
			this->vuDimLevels.clear();

			for(vector<size_t>::const_iterator 
				ivuDimLength = vuDimLengths.begin();
				ivuDimLength != vuDimLengths.end();
				ivuDimLength++)
			{
				size_t uDimLength = *ivuDimLength;
				size_t uDimLevel  = (size_t)ceil(log((double)uDimLength)/M_LN2) + 1;
				this->vuDimLevels.push_back(uDimLevel);
				size_t uCoefLength = 1 << (uDimLevel - 1);
				this->vuCoefLengths.push_back(uCoefLength);

				vector<size_t> vuSubLevel2Coef;
				vuSubLevel2Coef.resize(uCoefLength * uDimLevel);
				for(size_t p = 0, i = 0; i < uCoefLength; i++)
					for(size_t 
						l = 0, uW = 1 << uDimLevel, uAccumNrOfWins = 0, uNrOfWins = 1; 
						l < uDimLevel; 
						uAccumNrOfWins += uNrOfWins, uNrOfWins = (!l)?1:(uNrOfWins*2), l++, p++, uW >>= 1)
					{
						size_t uLevelWin = i / uW;
						vuSubLevel2Coef[p] = uAccumNrOfWins + uLevelWin;
					}

				this->vvuSubLevel2Coef.push_back(vuSubLevel2Coef);
			}

			int iWaveletDenomiator = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				iWaveletDenomiator *= 1 << (vuDimLevels[d] - 1);
			dWaveletDenomiator = sqrt((double)iWaveletDenomiator);

			size_t uThreshold = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				uThreshold *= this->vuCoefLengths[d];
			dWaveletThreshold = 1.0 / sqrt((double)uThreshold);

			_SetMaxLevels(this->vuDimLevels);
		}
		
		CSepDWTHeader()
		{
		}
	};
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
