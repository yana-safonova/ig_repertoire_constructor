#include "custom_vdj_hits_calculator.hpp"

using namespace vdj_labeler;

void CustomVDJHitsCalculator::AddHits(VDJHitsPtr vdj_hits, ImmuneGeneSegmentHitsPtr ig_gene_hits) {
    for(auto it = ig_gene_hits->cbegin(); it != ig_gene_hits->cend(); it++)
        vdj_hits->AddIgGeneAlignment(*it);
}

VDJHitsStoragePtr CustomVDJHitsCalculator::ComputeHits() const {
    VDJHitsStoragePtr vdj_hits_storage(new VDJHitsStorage(vj_alignment_info_));
    for(auto && vdj_hits : *vdj_hits_storage) {
        vdj_hits->SetDHits(d_hits_calculator_.ComputeDHits(vdj_hits->Read()));
    }
        /*
        cout << "==== Read: " << (*it)->name << ". # V alignments: " << vdj_hits_ptr->VHitsNumber() <<
                ", # D alignments: " << vdj_hits_ptr->DHitsNumber() << ", # J alignments: " <<
                vdj_hits_ptr->JHitsNumber() << endl;
        cout << "== V hits: " << endl;
        for(auto it = v_hits->cbegin(); it != v_hits->cend(); it++) {
            cout << **it << endl << endl;
        }
        cout << "== D hits" << endl;
        for(auto it = d_hits->cbegin(); it != d_hits->cend(); it++) {
            cout << **it << endl << endl;
        }
        cout << "== J hits" << endl;
        for(auto it = j_hits->cbegin(); it != j_hits->cend(); it++) {
            cout << **it << endl << endl;
        }
        */
    return vdj_hits_storage;
}
