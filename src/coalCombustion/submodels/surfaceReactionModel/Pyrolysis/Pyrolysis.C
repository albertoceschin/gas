#include "Pyrolysis.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Pyrolysis<CloudType>::
    Pyrolysis
(
 const dictionary& dict,
 CloudType& owner
 )
    :
        SurfaceReactionModel<CloudType>(dict, owner, typeName),
        Sb_(this->coeffDict().getScalar("Sb")),
        C1_(this->coeffDict().getScalar("C1")),
        C2_(this->coeffDict().getScalar("C2")),
        E_(this->coeffDict().getScalar("E")),
        C3_(this->coeffDict().getScalar("C3")),
        CsLocalId_(-1),
        C8H18LocalId_(-1),
        C10H22LocalId_(-1),
        O2GlobalId_(owner.composition().carrierId("O2")),
        H2GlobalId_(owner.composition().carrierId("H2")),
        CO2GlobalId_(owner.composition().carrierId("CO2")),
        C8H18GlobalId_(owner.composition().carrierId("C8H18")),
        C10H22GlobalId_(owner.composition().carrierId("C10H22")),
        WC_(0.0),
        WC8H18_(0.0),
        WC10H22_(0.0),
        WO2_(0.0),
        WH2_(0.0),
        HcCO2_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Determine Cs ids
    label idLiquid = owner.composition().idLiquid();
    C8H18LocalId_ = owner.composition().localId(idLiquid, "C8H18");
    C10H22LocalId_ = owner.composition().localId(idLiquid, "C10H22");

    // Set local copies of thermo properties
    WO2_ = owner.thermo().carrier().W(O2GlobalId_);
    WH2_ = owner.thermo().carrier().W(H2GlobalId_); 
    WC8H18_ = owner.thermo().carrier().W(C8H18GlobalId_); 
    WC10H22_ = owner.thermo().carrier().W(C10H22GlobalId_); 
    const scalar WCO2 = owner.thermo().carrier().W(CO2GlobalId_);
    WC_ = WCO2 - WO2_;

    HcCO2_ = owner.thermo().carrier().Hc(CO2GlobalId_);

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;
    const scalar YC10H22loc = owner.composition().Y0(idLiquid)[C10H22LocalId_];
    const scalar YLiquidTot1 = owner.composition().YMixture0()[idLiquid];
    Info<< "    C10H22(l): particle mass fraction = " << YC10H22loc*YLiquidTot1 << endl;
}


template<class CloudType>
Foam::Pyrolysis<CloudType>::
    Pyrolysis
(
 const Pyrolysis<CloudType>& srm
 )
    :
        SurfaceReactionModel<CloudType>(srm),
        Sb_(srm.Sb_),
        C1_(srm.C1_),
        C2_(srm.C2_),
        E_(srm.E_),
        C3_(srm.C3_),
        CsLocalId_(srm.CsLocalId_),
        C8H18LocalId_(srm.C8H18LocalId_),
        C10H22LocalId_(srm.C10H22LocalId_),
        O2GlobalId_(srm.O2GlobalId_),
        H2GlobalId_(srm.H2GlobalId_),
        CO2GlobalId_(srm.CO2GlobalId_),
        C8H18GlobalId_(srm.C8H18GlobalId_),
        C10H22GlobalId_(srm.C10H22GlobalId_),
        WC_(srm.WC_),
        WC8H18_(srm.WC8H18_),
        WC10H22_(srm.WC10H22_),
        WO2_(srm.WO2_),
        WH2_(srm.WH2_),
        HcCO2_(srm.HcCO2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::Pyrolysis<CloudType>::calculate
(
 const scalar dt,
 const scalar Re,
 const scalar nu,
 const label celli,
 const scalar d,
 const scalar T,
 const scalar Tc,
 const scalar pc,
 const scalar rhoc,
 const scalar mass,
 const scalarField& YGas,
 const scalarField& YLiquid,
 const scalarField& YSolid,
 const scalarField& YMixture,
 const scalar N,
 scalarField& dMassGas,
 scalarField& dMassLiquid,
 scalarField& dMassSolid,
 scalarField& dMassSRCarrier
 ) const
{

    scalar dh (0.);

    const label idSolid = CloudType::parcelType::SLD;
    const label idLiquid = CloudType::parcelType::LIQ;

    const scalar fCombCs = YMixture[idSolid]*YSolid[CsLocalId_];
    const scalar fCombC10H22 = YMixture[idLiquid]*YLiquid[C10H22LocalId_];

    if (fCombCs > SMALL)
    {
        const SLGThermo& thermo = this->owner().thermo();
        const scalar YO2 = thermo.carrier().Y(O2GlobalId_)[celli];
        const scalar D0 = C1_/d*pow(0.5*(T + Tc), 0.75);
        const scalar Rk = C2_*exp(-E_/(RR*Tc));
        const scalar Ap = constant::mathematical::pi*sqr(d);
        scalar dmC = Ap*rhoc*RR*Tc*YO2/WO2_*D0*Rk/(D0 + Rk)*dt;
        //Info<<mass*fCombCs << "..."<< dmC<<endl;
        dmC = min(mass*fCombCs, dmC);
        const scalar dOmega = dmC/WC_;
        //Info<<dOmega<<"dOmega"<<endl;
        const scalar dmO2 = dOmega*Sb_*WO2_;
        const scalar dmCO2 = dOmega*(WC_ + Sb_*WO2_);
        dMassSolid[CsLocalId_] += dOmega*WC_;
        dMassSRCarrier[O2GlobalId_] -= dmO2;
        dMassSRCarrier[CO2GlobalId_] += dmCO2;
        const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);
        // Heat of reaction [J]
        dh += dmC*HsC - dmCO2*HcCO2_;
    }

    if (fCombC10H22 > SMALL)
    {
        scalar dmC10H22 =  C3_ *dt;
        dmC10H22 = min(mass*fCombC10H22, dmC10H22);
        const scalar dOmega = dmC10H22/WC10H22_;
        const scalar dmH2 = dOmega*2.*WH2_;
        dMassSolid[CsLocalId_] -= dOmega*WC_;
        dMassLiquid[C10H22LocalId_] += dOmega*WC10H22_;
        dMassLiquid[C8H18LocalId_] -= dOmega*WC8H18_;
        dMassSRCarrier[H2GlobalId_] += dmH2;
    }

    return dh;
}


// ************************************************************************* //
