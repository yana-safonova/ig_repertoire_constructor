#include "custom_vdj_hits_calculator.hpp"

using namespace std;

void CustomVDJHitsCalculator::AddHits(VDJHitsPtr vdj_hits, IgGeneSegmentHitsPtr ig_gene_hits) {
    for(auto it = ig_gene_hits->cbegin(); it != ig_gene_hits->cend(); it++)
        vdj_hits->AddIgGeneAlignment(*it);
}

VDJHitsStoragePtr CustomVDJHitsCalculator::ComputeHits() {
    VDJHitsStoragePtr vdj_hits_storage(new VDJHitsStorage());
    for(auto it = read_archive_.cbegin(); it != read_archive_.cend(); it++) {
        if(!vj_alignment_info_.ContainsRead((*it)->id))
            continue;
        IgGeneSegmentHitsPtr v_hits = v_hits_calculator_.ComputeHits(*it);
        IgGeneSegmentHitsPtr j_hits = j_hits_calculator_.ComputeHits(*it);
        IgGeneSegmentHitsPtr d_hits = d_hits_calculator_.ComputeHits(*it);
        VDJHitsPtr vdj_hits_ptr(new VDJHits(*it));
        AddHits(vdj_hits_ptr, v_hits);
        AddHits(vdj_hits_ptr, j_hits);
        AddHits(vdj_hits_ptr, d_hits);
        vdj_hits_storage->AddVDJHits(vdj_hits_ptr);
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
    }
    return vdj_hits_storage;
}
