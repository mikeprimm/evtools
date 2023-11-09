import { TestBed } from '@angular/core/testing';

import { ExoplanetFetchService } from './exoplanet-fetch.service';

describe('ExoplanetFetchService', () => {
  let service: ExoplanetFetchService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(ExoplanetFetchService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
