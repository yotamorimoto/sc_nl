// NLUGens
// http://yota.tehis.net/

#include "SC_PlugIn.h"
#include <limits>

#define LATTICE 10
#define GENEBIT 32

static InterfaceTable *ft;

struct Logist : public Unit {
	double x;
	float counter;
};
struct Nagumo : public Unit {
	double u, v;
};
struct FIS : public Unit {
};
struct CML : public Unit {
	double x[LATTICE];
	float counter;
};
struct GCM : public Unit {
	double x[LATTICE];
	float counter;
};
struct HCM : public Unit {
	uint32 x[GENEBIT];
	float counter;
};
struct TLogist : public Logist {
	double trig;
};

static void Logist_next(Logist *unit, int inNumSamples);
static void Logist_Ctor(Logist *unit);
static void Nagumo_next(Nagumo *unit, int inNumSamples);
static void Nagumo_Ctor(Nagumo *unit);
static void FIS_next(FIS *unit, int inNumSamples);
static void FIS_Ctor(FIS *unit);
static void CML_next(CML *unit, int inNumSamples);
static void CML_Ctor(CML *unit);
static void GCM_next(GCM *unit, int inNumSamples);
static void GCM_Ctor(GCM *unit);
static void HCM_next(HCM *unit, int inNumSamples);
static void HCM_Ctor(HCM *unit);
static void TLogist_next(TLogist *unit, int inNumSamples);
static void TLogist_Ctor(TLogist *unit);

inline double logist(double r, double x);
inline double logist(double r, double x)
{
	return 1.l - r * x * x;
}

void Logist_next(Logist *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);

	double x = unit->x;
	float counter = unit->counter;

	float spc;
	if(freq < SAMPLERATE)
		spc = SAMPLERATE / sc_max(freq, 0.001f);
	else spc = 1.f;

	LOOP(inNumSamples,
		 if(counter >= spc){
			 counter -= spc;
			 x = logist(r, x);
		 }
		 counter++;
		 ZXP(out) = x;
		 )
	unit->x = x;
	unit->counter = counter;
}

void Logist_Ctor(Logist *unit)
{
	SETCALC(Logist_next);
	unit->x = IN0(2);
	unit->counter = 0.f;
	Logist_next(unit, 1);
}


void Nagumo_next(Nagumo *unit, int inNumSamples)
{
	float *out = ZOUT(0);

	double uh = ZIN0(0);
	double vh = ZIN0(1);
	float *pulse = ZIN(2);

	double u = unit->u;
	double v = unit->v;

	LOOP(inNumSamples,
		float zPulse = ZXP(pulse);
		u += uh * (10.l * (- v + u - 0.3333333l*u*u*u + zPulse));
		v += vh * (u - 0.8l * v + 0.7l);
		ZXP(out) = u * 0.3;
	)
	unit->u = u;
	unit->v = v;
}

void Nagumo_Ctor(Nagumo *unit)
{
	SETCALC(Nagumo_next);
	unit->u = 0.1l;
	unit->v = 0.l;
	Nagumo_next(unit, 1);
}

void FIS_next(FIS *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float *r = ZIN(0);
	float *x = ZIN(1);
	int n = ZIN0(2);

	LOOP(inNumSamples,
		double zx = ZXP(x);
		double zr = ZXP(r);
		for(int i=0; i<n; i++)
			zx = sin(zr * zx);
		ZXP(out) = zx;
	)
}

void FIS_Ctor(FIS *unit)
{
	SETCALC(FIS_next);
	FIS_next(unit, 1);
}

void CML_next(CML *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	double x[LATTICE];

	memcpy(x, unit->x, sizeof(unit->x));
	float counter = unit->counter;

	float spc;
	if (freq < SAMPLERATE)
		spc = SAMPLERATE / sc_max(freq, 0.001f);
	else spc = 1.f;

	LOOP(inNumSamples,
		if (counter >= spc) {
			counter -= spc;
			for (int i=0; i<LATTICE; i++) {
				// one-way coupling sonically more interesting
				// also avoids negative modulo
				x[i] = (1.l - g) * logist(r, x[i]) + g * logist(r, x[(i+1)%LATTICE]);
			}
		}
		counter++;
		ZXP(out) = x[5];
	)
	unit->counter = counter;
	memcpy(unit->x, x, sizeof(x));
}

void CML_Ctor(CML *unit)
{
	RGET
	SETCALC(CML_next);
	for (int i=0; i<LATTICE; i++) unit->x[i] = frand(s1,s2,s3);
	unit->counter = 0.f;
	CML_next(unit, 1);
	RPUT
}


void GCM_next(GCM *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	double x[LATTICE];
	double reciprocal = 1.l / LATTICE;
	float counter = unit->counter;
	memcpy(x, unit->x, sizeof(unit->x));

	float spc;
	if (freq < SAMPLERATE)
		spc = SAMPLERATE / sc_max(freq, 0.001f);
	else spc = 1.f;

	// once per block is more interesting
	double sum = 0.l;
	for (int i=0; i<LATTICE; i++) sum += logist(r, x[i]);

	LOOP(inNumSamples,
		if (counter >= spc) {
			counter -= spc;
			for (int i=0; i<LATTICE; i++) {
				x[i] = (1.l - g) * logist(r, x[i]) + g * reciprocal * sum;
			}
		}
		counter++;
		ZXP(out) = x[5];
	)
	unit->counter = counter;
	memcpy(unit->x, x, sizeof(x));
}

void GCM_Ctor(GCM *unit)
{
	RGET
	SETCALC(GCM_next);
	for (int i=0; i<LATTICE; i++) unit->x[i] = frand(s1,s2,s3);
	unit->counter = 0.f;
	GCM_next(unit, 1);
	RPUT
}


inline unsigned flip(unsigned x, unsigned bit);
inline unsigned flip(unsigned x, unsigned bit)
{
	return x ^ (1UL << bit);
}
// uint32
uint32 halfmax = std::numeric_limits<uint32>::max() / 2;
inline double l2f(uint32 in);
inline double l2f(uint32 in)
{
	return (double)in / halfmax - 1.l;
}
inline uint32 f2l(double in);
inline uint32 f2l(double in)
{
	return (uint32)(in * halfmax + halfmax);
}
void HCM_next(HCM *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	uint32 x[GENEBIT];
	double reciprocal = 1.l / GENEBIT;

	memcpy(x, unit->x, sizeof(unit->x));
	float counter = unit->counter;

	float spc;
	if (freq < SAMPLERATE)
		spc = SAMPLERATE / sc_max(freq, 0.001f);
	else spc = 1.f;

	// once per control is more interesting
	// provided we have sufficient bit space like uint32
	// otherwise it will repeat (which maybe somewhat useful?)
	double sum = 0.l;
	for (int i=0; i<GENEBIT; i++)
		for (int j=0; j<GENEBIT; j++) sum += logist(r, l2f(flip(x[j], i)));

	LOOP(inNumSamples,
		 if (counter >= spc) {
			 counter -= spc;
			 for (int i=0; i<GENEBIT; i++) {
				 // uint16 is enough if we renew the sum everytime?
				 // for (int j=0; j<GENEBIT; j++) sum += logist(r, l2f(flip(x[j], i)));
				 double tmp = (1.l - g) * logist(r, l2f(x[i])) + g * reciprocal * sum;
				 x[i] = f2l(tmp);
			 }
		 }
		 counter++;
		 ZXP(out) = l2f(x[15]);
	)
	memcpy(unit->x, x, sizeof(x));
	unit->counter = counter;
}

void HCM_Ctor(HCM *unit)
{
	SETCALC(HCM_next);
	for (int i=0; i<GENEBIT; i++) unit->x[i] = i%2;
	unit->counter = 0.f;
	HCM_next(unit, 1);
}


void TLogist_next(TLogist *unit, int inNumSamples)
{
	float trig = ZIN0(2);
	if (trig > 0.f && unit->trig <= 0.f) {
		double r = ZIN0(0);
		ZOUT0(0) = unit->x = 1.f - r * unit->x * unit->x;
	} else {
		ZOUT0(0) = unit->x;
	}
	unit->trig = trig;
}

void TLogist_Ctor(TLogist *unit)
{
	double r = ZIN0(0);
	ZOUT0(0) = unit->x = ZIN0(1);
	SETCALC(TLogist_next);
	unit->trig = ZIN0(2);
}

PluginLoad(NL)
{
	ft = inTable;
	DefineSimpleUnit(Logist);
	DefineSimpleUnit(Nagumo);
	DefineSimpleUnit(FIS);
	DefineSimpleUnit(CML);
	DefineSimpleUnit(GCM);
	DefineSimpleUnit(HCM);
	DefineSimpleUnit(TLogist);
}
