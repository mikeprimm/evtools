import { ComponentFixture, TestBed } from '@angular/core/testing';

import { ExoplanetsComponent } from './exoplanets.component';

describe('ExoplanetsComponent', () => {
  let component: ExoplanetsComponent;
  let fixture: ComponentFixture<ExoplanetsComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ ExoplanetsComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(ExoplanetsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
